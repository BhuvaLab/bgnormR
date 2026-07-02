#' Fit a Gaussian Mixture Model to intensity data
#'
#' Wraps \code{mclust::Mclust} with the variable-variance (\code{"V"}) model
#' and returns a tidy list of parameters and posterior probabilities for the
#' fitting subset.  Use \code{\link{.gmm_posteriors}} to obtain posteriors for
#' the full input vector.
#'
#' @param x Numeric vector of (log-transformed) intensities.
#' @param n_components Integer, number of GMM components (2 or 3).
#' @param sample_prop Numeric in (0, 1]; proportion of non-zero observations to
#'   use for fitting.  Values < 1 draw a random subset without replacement.
#'   Default \code{1} (use all non-zero observations).
#' @param ... Extra arguments forwarded to \code{mclust::Mclust}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{means}{Component means (length \code{n_components}).}
#'     \item{sds}{Component standard deviations.}
#'     \item{props}{Mixing proportions.}
#'     \item{posteriors}{Matrix (n_fit x G) of posterior probabilities for the
#'       fitting subset only.  Use \code{.gmm_posteriors()} for the full data.}
#'   }
#' @import mclust
#' @keywords internal
.fit_gmm <- function(x, n_components = 3L, sample_prop = 1,
                     .sample_threshold = 1e5, ...) {
  stopifnot(n_components %in% 2:3)
  stopifnot(is.numeric(sample_prop), length(sample_prop) == 1L,
            sample_prop > 0, sample_prop <= 1)

  # Remove zeros (unexposed / masked pixels) and non-finite values
  x <- x[is.finite(x) & x != 0]
  n_nonzero <- length(x)

  # Subsample only when non-zero pixel count exceeds the threshold; otherwise
  # use all data and inform the user that sampling is not needed.
  if (sample_prop < 1) {
    if (n_nonzero > .sample_threshold) {
      n_keep <- max(n_components * 10L, as.integer(ceiling(sample_prop * n_nonzero)))
      x <- x[sample.int(n_nonzero, min(n_keep, n_nonzero), replace = FALSE)]
    } else {
      message("Sampling not required: using all ", n_nonzero,
              " non-zero pixels for GMM fitting.")
    }
  }

  if (length(x) < n_components * 10L) {
    stop("Too few finite non-zero observations to fit a ", n_components,
         "-component GMM (need >= ", n_components * 10L, ").")
  }

  fit <- mclust::Mclust(x, G = n_components, modelNames = "V",
                        verbose = FALSE, ...)

  if (is.null(fit)) {
    stop("GMM fitting failed (mclust returned NULL). ",
         "Check that the data has sufficient variation.")
  }

  # Extract variance: for modelNames="V" the variance is per-component
  var_vec <- if (is.matrix(fit$parameters$variance$sigmasq)) {
    diag(fit$parameters$variance$sigmasq)
  } else {
    fit$parameters$variance$sigmasq
  }

  list(
    means      = fit$parameters$mean,
    sds        = sqrt(var_vec),
    props      = fit$parameters$pro,
    posteriors = fit$z   # for the fitting subset only
  )
}

#' Compute GMM posterior probabilities for an arbitrary input vector
#'
#' Given fitted GMM parameters (from \code{\link{.fit_gmm}}), computes the
#' posterior probability of each observation belonging to each component.
#' This should be called on the full data vector after fitting on a subsample.
#'
#' @param x     Numeric vector; the full set of observations.
#' @param means Component means.
#' @param props Mixing proportions.
#' @param sds   Component standard deviations.
#'
#' @return A matrix of dimension \code{length(x) x G} with rows summing to 1.
#' @keywords internal
.gmm_posteriors <- function(x, means, props, sds) {
  G <- length(means)
  # Weighted density per component
  dens <- vapply(seq_len(G), function(k) {
    props[k] * stats::dnorm(x, means[k], sds[k])
  }, numeric(length(x)))
  # Normalise rows to sum to 1
  row_tots <- rowSums(dens)
  row_tots[row_tots == 0] <- .Machine$double.eps
  dens / row_tots
}

#' Classify observations using threshold-corrected GMM assignment
#'
#' Starts from the posterior-probability argmax, then applies two correction
#' passes for the 3-component model so that the distribution lower tails never
#' overwhelm the hard evidence of a low pixel value:
#' \enumerate{
#'   \item If argmax = 3 (Signal) but \code{x < means[2]} -> demote to 2.
#'   \item If class is now 2 (Non-specific) but \code{x < means[1]} -> demote to 1.
#' }
#' For the 2-component model the argmax is returned unchanged.
#'
#' @param x         Numeric vector of (log-transformed) intensities, same
#'   length as \code{nrow(posteriors)}.
#' @param means     Resolved component means (1 = background, 2 = non-specific,
#'   3 = signal).
#' @param posteriors Matrix (n x G) of posterior probabilities.
#'
#' @return Integer vector of component labels (values in \code{1:G}).
#' @keywords internal
.classify_pixels <- function(x, means, posteriors) {
  G   <- length(means)
  cls <- max.col(posteriors, ties.method = "first")
  if (G == 3L) {
    cls[cls == 3L & x < means[2L]] <- 2L   # Signal -> Non-specific if below mean_2
    cls[cls == 2L & x < means[1L]] <- 1L   # Non-specific -> Background if below mean_1
  }
  cls
}

#' Compute the Non-specific / Signal class boundary threshold
#'
#' Given a GMM-based per-pixel classification (from \code{.classify_pixels})
#' and the corresponding background-adjusted intensities, computes the maximum
#' adjusted value among pixels classified as Non-specific (class 2).  This
#' single scalar is sufficient to recover the full 3-class assignment:
#' \itemize{
#'   \item Class 1 (Background):   \code{adjusted == 0}
#'   \item Class 2 (Non-specific): \code{0 < adjusted <= threshold}
#'   \item Class 3 (Signal):       \code{adjusted > threshold}
#' }
#' The class-1 boundary is always 0 (enforced by the \code{pmax} in
#' \code{.deconvolve}) and is not stored.
#'
#' @param x_adj Numeric vector of background-adjusted intensities (length n,
#'   all values \code{>= 0} after \code{pmax}).
#' @param cls   Integer vector of class labels produced by
#'   \code{.classify_pixels} (values in \code{1:3}).
#'
#' @return A single numeric scalar: the maximum adjusted intensity of
#'   class-2 pixels, or \code{0} if no pixels are classified as Non-specific
#'   (in which case every non-zero pixel is Signal by convention).
#' @keywords internal
.compute_threshold <- function(x_adj, cls) {
  v <- x_adj[cls == 2L]
  if (length(v) > 0L) max(v) else 0
}

#' Precompute a histogram of log-transformed intensities for storage
#'
#' Used during \code{bgnorm_pixels} and \code{bgnorm_cells} fitting to store a
#' compact histogram of the raw data so that \code{plot_distributions} does not
#' need the original image or assay.
#'
#' @param x_log Numeric vector of log-transformed intensities (may contain
#'   zeros and non-finite values; these are excluded before binning).
#' @param bins  Integer number of histogram bins (default 100).
#'
#' @return A list with elements \code{breaks} (length bins + 1) and
#'   \code{density} (length bins), as returned by \code{hist(..., plot = FALSE)}.
#' @keywords internal
.compute_histogram <- function(x_log, bins = 100L) {
  ok   <- is.finite(x_log) & x_log != 0
  x_ok <- x_log[ok]
  rng  <- range(x_ok)
  brks <- seq(rng[1L], rng[2L], length.out = as.integer(bins) + 1L)
  h    <- hist(x_ok, breaks = brks, plot = FALSE)
  list(breaks = h$breaks, density = h$density)
}

#' Fit G=2 and G=3 GMMs on the same pixel sample and compare BIC
#'
#' Filters zeros, optionally subsamples, then calls \code{mclust::Mclust} with
#' \code{G = 2:3} so that both models are fitted on an identical dataset and
#' their BIC values are directly comparable.  mclust selects the higher-BIC
#' model; \code{no_signal} is \code{TRUE} when G=2 is selected.
#'
#' @param x         Numeric vector of log-transformed intensities.
#' @param sample_prop Numeric in (0, 1]; subsampling fraction (default 1).
#' @param .sample_threshold Integer; only subsample when \code{n_nonzero}
#'   exceeds this value (default 1e5).
#' @param ...       Extra arguments forwarded to \code{mclust::Mclust}.
#'
#' @return A list with elements \code{means}, \code{sds}, \code{props},
#'   \code{posteriors} (fitting subset, G components), \code{bic} (named
#'   \code{c(G2 = ..., G3 = ...)}, each divided by the number of points
#'   used to fit so values are point-specific and comparable across channels),
#'   and \code{no_signal} (logical).
#' @keywords internal
.fit_gmm_bic <- function(x, sample_prop = 1, .sample_threshold = 1e5, ...) {
  x <- x[is.finite(x) & x != 0]
  n_nonzero <- length(x)

  if (sample_prop < 1) {
    if (n_nonzero > .sample_threshold) {
      n_keep <- max(30L, as.integer(ceiling(sample_prop * n_nonzero)))
      x <- x[sample.int(n_nonzero, min(n_keep, n_nonzero), replace = FALSE)]
    } else {
      message("Sampling not required: using all ", n_nonzero,
              " non-zero pixels for GMM fitting.")
    }
  }

  n_fit <- length(x)
  if (n_fit < 30L)
    stop("Too few finite non-zero observations to fit GMM (need >= 30).")

  fit <- mclust::Mclust(x, G = 2:3, modelNames = "V", verbose = FALSE, ...)

  if (is.null(fit))
    stop("GMM fitting failed (mclust returned NULL). ",
         "Check that the data has sufficient variation.")

  bic_g2 <- fit$BIC["2", "V"] / n_fit
  bic_g3 <- fit$BIC["3", "V"] / n_fit

  var_vec <- if (is.matrix(fit$parameters$variance$sigmasq)) {
    diag(fit$parameters$variance$sigmasq)
  } else {
    fit$parameters$variance$sigmasq
  }

  list(
    means      = fit$parameters$mean,
    sds        = sqrt(var_vec),
    props      = fit$parameters$pro,
    posteriors = fit$z,
    bic        = c(G2 = bic_g2, G3 = bic_g3),
    no_signal  = (fit$G == 2L)
  )
}

#' Resolve GMM component identities
#'
#' Orders components as background (1), non-specific/autofluorescence (2),
#' and biological signal (3).  When the two highest-mean components are
#' close, the one with the heavier right tail (higher 75th percentile) is
#' designated signal, following Kharbanda et al. 2025.
#'
#' @param means Numeric vector of component means.
#' @param sds   Numeric vector of component standard deviations.
#' @param props Numeric vector of mixing proportions.
#' @param posteriors Matrix of posterior probabilities (n x G).
#'
#' @return A list with \code{means}, \code{sds}, \code{props}, and
#'   \code{posteriors} reordered so component 1 = background, 2 = nonspecific,
#'   3 = signal.
#' @keywords internal
.resolve_components <- function(means, sds, props, posteriors) {
  G <- length(means)
  ord <- order(means)          # ascending by mean: 1=background, G=signal candidate

  if (G == 3L) {
    i2 <- ord[2L]; i3 <- ord[3L]
    # Compare right tails: 75th percentile
    q75_2 <- stats::qnorm(0.75, means[i2], sds[i2])
    q75_3 <- stats::qnorm(0.75, means[i3], sds[i3])
    # If candidate-2 has heavier right tail than candidate-3, swap
    if (q75_2 > q75_3) {
      ord[2L:3L] <- ord[c(3L, 2L)]
    }
  }

  list(
    means      = means[ord],
    sds        = sds[ord],
    props      = props[ord],
    posteriors = posteriors[, ord, drop = FALSE]
  )
}
