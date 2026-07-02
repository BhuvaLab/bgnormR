#' Jensen-Shannon Divergence between two Gaussian distributions
#'
#' Computes JSD numerically on a dense grid, used as a signal-to-noise QC
#' metric between the non-specific binding component (2) and the biological
#' signal component (3) of the bgnorm model.  Higher values indicate
#' better separation and staining quality.
#'
#' @param mu1,mu2 Means of the two Gaussian distributions.
#' @param sd1,sd2 Standard deviations of the two Gaussian distributions.
#' @param n_grid  Number of evaluation points (default 2000).
#'
#' @return Scalar JSD value in [0, 1] (bits, log2 base).
#' @export
#' @examples
#' # Well-separated distributions -> high JSD
#' jsd_gaussians(mu1 = 0, sd1 = 0.5, mu2 = 3, sd2 = 0.8)
#' # Identical distributions -> JSD = 0
#' jsd_gaussians(mu1 = 1, sd1 = 1, mu2 = 1, sd2 = 1)
jsd_gaussians <- function(mu1, sd1, mu2, sd2, n_grid = 2000L) {
  stopifnot(is.numeric(c(mu1, sd1, mu2, sd2)), sd1 > 0, sd2 > 0)

  lo <- min(mu1 - 5 * sd1, mu2 - 5 * sd2)
  hi <- max(mu1 + 5 * sd1, mu2 + 5 * sd2)
  x  <- seq(lo, hi, length.out = n_grid)
  dx <- x[2L] - x[1L]

  p1 <- stats::dnorm(x, mu1, sd1)
  p2 <- stats::dnorm(x, mu2, sd2)
  m  <- 0.5 * (p1 + p2)

  eps <- .Machine$double.eps
  kl1 <- sum(p1 * log2((p1 + eps) / (m + eps)) * dx)
  kl2 <- sum(p2 * log2((p2 + eps) / (m + eps)) * dx)

  0.5 * (kl1 + kl2)
}

#' Compute the JSD quality control metric from bgnorm model parameters
#'
#' Extracts the JSD between components 2 (non-specific binding) and 3
#' (signal) from the fitted model parameters stored in a
#' \code{\link{BgnormResult}} object.
#'
#' @param result A \code{BgnormResult} object returned by
#'   \code{\link{bgnorm_pixels}} or \code{\link{bgnorm_cells}}.
#'
#' @return Named numeric scalar \code{jsd}.
#' @export
#' @examples
#' set.seed(1)
#' x   <- c(exp(rnorm(300, log(300), 0.5)), exp(rnorm(100, log(2000), 0.9)))
#' res <- bgnorm_cells(x)
#' jsd_qc(res)
jsd_qc <- function(result) {
  stopifnot(inherits(result, "BgnormResult"))
  p <- result$parameters
  G <- length(p$means)

  if (G == 3L) {
    jsd_gaussians(p$means[2L], p$sds[2L], p$means[3L], p$sds[3L])
  } else {
    # 2-component (cell level)
    jsd_gaussians(p$means[1L], p$sds[1L], p$means[2L], p$sds[2L])
  }
}

#' Summarise QC across multiple bgnorm results
#'
#' @param results A named list of \code{BgnormResult} objects (e.g., one per
#'   marker), or a \code{QPTIFFImage} returned by \code{\link{bgnorm_pixels}}
#'   (bgnorm results are extracted automatically).
#'
#' @return A \code{data.frame} with columns \code{name}, \code{jsd}, and
#'   \code{prop_signal} (proportion of pixels/cells in the signal component).
#' @export
#' @examples
#' set.seed(2)
#' mk1 <- bgnorm_cells(c(exp(rnorm(300,log(300),0.5)), exp(rnorm(100,log(2000),0.9))))
#' mk2 <- bgnorm_cells(c(exp(rnorm(400,log(250),0.4)), exp(rnorm( 80,log(1800),0.8))))
#' qc_summary(list(CD3 = mk1, CD8 = mk2))
qc_summary <- function(results) {
  if (inherits(results, "QPTIFFImage")) {
    results <- bgnorm_results(results)
    if (is.null(results))
      stop("'results' QPTIFFImage has no bgnorm results. ",
           "Pass the QPTIFFImage returned by bgnorm_pixels().")
  }
  stopifnot(is.list(results),
            all(vapply(results, inherits, logical(1L), "BgnormResult")))

  nms <- names(results)
  if (is.null(nms)) nms <- paste0("item", seq_along(results))

  rows <- lapply(seq_along(results), function(i) {
    r  <- results[[i]]
    G  <- length(r$parameters$means)
    jsd_val  <- r$jsd
    prop_sig <- if (isTRUE(r$no_signal)) NA_real_ else r$parameters$props[G]
    data.frame(name = nms[i], jsd = jsd_val, prop_signal = prop_sig,
               stringsAsFactors = FALSE)
  })

  do.call(rbind, rows)
}
