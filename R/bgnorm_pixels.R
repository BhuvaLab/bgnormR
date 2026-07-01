#' Pixel-level background normalisation (bgnorm) for a QPTIFFImage
#'
#' Fits a three-component Gaussian Mixture Model (GMM) to the log2-transformed
#' pixel intensities for each channel in a \code{\link{QPTIFFImage}} and
#' deconvolves the biological signal component from background sources.  All
#' channels are processed independently; use \code{channels} to restrict
#' processing to a subset.  Optionally runs channels in parallel via
#' \pkg{BiocParallel}.
#'
#' The preprocessing, model, and deconvolution steps follow Kharbanda et al.
#' 2025:
#' \enumerate{
#'   \item Log2-transform: \eqn{I_{\log} = \log_2(I/c + 1)}.
#'   \item Fit a 3-component GMM via EM (mclust), filtering zeros and
#'         optionally subsampling for speed.
#'   \item Resolve component identities (background, non-specific, signal).
#'   \item Apply the deconvolution formula to recover \eqn{U_3}.
#'   \item Optionally normalise by the 75th-percentile of the signal
#'         component (bgnormQ).
#' }
#'
#' @param img         A \code{\link{QPTIFFImage}} \code{[H x W x C]}.  Use
#'   \code{\link{as.QPTIFFImage}} to promote a plain 3-D array or 2-D matrix.
#' @param channels    Character or integer vector of channels to process, or
#'   \code{NULL} (default) to process all channels.
#' @param cofactor    Positive numeric cofactor for the log2 transform.
#'   Default 150 is appropriate for 16-bit Akoya PhenoCycler-Fusion images.
#' @param quantile_norm Logical; apply model-based quantile normalisation
#'   (bgnormQ)?  Default \code{FALSE}.
#' @param quantile    Quantile of the signal component used for normalisation
#'   when \code{quantile_norm = TRUE}.  Default 0.75.
#' @param sample_prop Numeric in (0, 1]; proportion of non-zero pixels used to
#'   fit the GMM.  Values < 1 draw a random subset, accelerating fitting on
#'   large images.  Default \code{1} (all non-zero pixels).
#' @param BPPARAM     A \code{\link[BiocParallel]{BiocParallelParam}} instance.
#'   Default \code{BiocParallel::SerialParam()}.
#' @param ...         Additional arguments forwarded to the internal GMM
#'   fitter.
#'
#' @return A \code{QPTIFFImage} \code{[H x W x C]} of background-adjusted
#'   log-intensities.  The per-channel \code{\link{BgnormResult}} objects
#'   (GMM parameters, posteriors, JSD) are stored as an attribute; retrieve
#'   them with \code{\link{bgnorm_results}(result)}.
#'
#' @seealso \code{\link{bgnorm_cells}}, \code{\link{bgnorm_results}},
#'   \code{\link{jsd_qc}}, \code{\link{as.QPTIFFImage}}
#' @importFrom BiocParallel bplapply SerialParam
#' @export
#' @examples
#' set.seed(42)
#' px  <- c(exp(rnorm(500, log(50),   0.4)),
#'          exp(rnorm(400, log(300),  0.6)),
#'          exp(rnorm(100, log(2000), 1.0)))
#' arr <- array(c(px, px), dim = c(length(px), 1L, 2L),
#'              dimnames = list(NULL, NULL, c("DAPI", "CD3")))
#' img <- as.QPTIFFImage(arr)
#' res <- bgnorm_pixels(img)
#' names(res)                        # channel names of the adjusted image
#' bgnorm_results(res)[["DAPI"]]     # BgnormResult for the DAPI channel
#' # Process only selected channels
#' res2 <- bgnorm_pixels(img, channels = "DAPI")
#'
#' \donttest{
#' # Real tissue image shipped with the package
#' path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' res  <- bgnorm_pixels(img, sample_prop = 0.1)
#' bgnorm_results(res)[["PanCK"]]    # per-channel model summary
#' }
bgnorm_pixels <- function(img, channels = NULL, cofactor = 150,
                           quantile_norm = FALSE, quantile = 0.75,
                           sample_prop = 0.1,
                           BPPARAM = BiocParallel::SerialParam(), ...) {
  if (!inherits(img, "QPTIFFImage"))
    stop("'img' must be a QPTIFFImage. ",
         "Use as.QPTIFFImage() to convert a 3-D array or 2-D matrix.")

  # Channel selection
  all_chs <- names(img)
  if (!is.null(channels)) {
    if (is.numeric(channels)) channels <- all_chs[as.integer(channels)]
    miss <- setdiff(channels, all_chs)
    if (length(miss))
      stop("Channels not found: ", paste(miss, collapse = ", "),
           "\nAvailable: ", paste(all_chs, collapse = ", "))
    img <- img[, , channels, drop = FALSE]
  }

  d        <- dim(img)
  chs      <- names(img)
  img_meta <- attr(img, "metadata") %||% list()

  mat <- matrix(as.array(img), nrow = d[1L] * d[2L], ncol = d[3L])
  colnames(mat) <- chs

  results <- .bgnorm_columns(
    x       = mat,
    fn      = function(col) .bgnorm_pixel_col(
      col,
      cofactor      = cofactor,
      quantile_norm = quantile_norm,
      quantile      = quantile,
      sample_prop   = sample_prop, ...
    ),
    BPPARAM = BPPARAM
  )

  # Reshape adjusted values into a QPTIFFImage [H x W x C].
  # The per-channel BgnormResult objects are stored as the bgnorm_results slot;
  # the QPTIFFImage itself holds the adjusted intensities.
  adj          <- attr(results, "adjusted_matrix")
  arr          <- array(adj,
                        dim      = c(d[1L], d[2L], ncol(adj)),
                        dimnames = list(NULL, NULL, colnames(adj)))
  results_list <- results
  attr(results_list, "adjusted_matrix") <- NULL

  .new_QPTIFFImage(arr, img_meta, bgnorm_results = results_list)
}

# ---- private: per-channel GMM algorithm ------------------------------------

.bgnorm_pixel_col <- function(x, cofactor = 150, quantile_norm = FALSE,
                               quantile = 0.75, sample_prop = 0.1, ...) {
  if (length(x) < 30L)
    stop("Too few observations to fit a GMM (need >= 30 pixels per channel).")

  x_log <- log_transform(x, cofactor = cofactor)

  # Fit G=2 and G=3 on the same pixel sample and compare BIC.
  # mclust selects the model with higher BIC; no_signal = TRUE when G=2 wins.
  gmm <- .fit_gmm_bic(x_log, sample_prop = sample_prop, ...)

  if (gmm$no_signal) {
    post2 <- .gmm_posteriors(x_log, gmm$means, gmm$props, gmm$sds)
    x_adj <- numeric(length(x_log))
    return(list(
      adjusted = x_adj,
      result   = .new_BgnormResult(
        parameters    = list(means = gmm$means, sds = gmm$sds, props = gmm$props),
        posteriors    = post2,
        threshold     = NULL,
        histogram     = .compute_histogram(x_log),
        jsd           = NA_real_,
        bic           = gmm$bic,
        level         = "pixel",
        quantile_norm = FALSE,
        no_signal     = TRUE
      )
    ))
  }

  full_post <- .gmm_posteriors(x_log, gmm$means, gmm$props, gmm$sds)
  params    <- .resolve_components(gmm$means, gmm$sds, gmm$props, full_post)
  x_adj     <- .deconvolve(x_log, params$means, params$sds, params$posteriors)

  if (quantile_norm) {
    q_factor <- .qnorm_factor(params$means, params$sds, params$props,
                               quantile = quantile)
    x_adj <- x_adj / q_factor
  }

  jsd <- jsd_gaussians(params$means[2L], params$sds[2L],
                        params$means[3L], params$sds[3L])

  cls <- .classify_pixels(x_log, params$means, params$posteriors)
  list(
    adjusted = x_adj,
    result   = .new_BgnormResult(
      parameters    = list(means = params$means, sds = params$sds,
                           props = params$props),
      posteriors    = params$posteriors,
      threshold     = .compute_threshold(x_adj, cls),
      histogram     = .compute_histogram(x_log),
      jsd           = jsd,
      bic           = gmm$bic,
      level         = "pixel",
      quantile_norm = quantile_norm
    )
  )
}

# ---- private helper: parallelise bgnorm over matrix columns ----------------

#' @keywords internal
.bgnorm_columns <- function(x, fn, BPPARAM) {
  markers <- colnames(x)
  if (is.null(markers)) markers <- paste0("marker", seq_len(ncol(x)))

  results <- BiocParallel::bplapply(
    seq_len(ncol(x)),
    function(j) tryCatch(fn(x[, j]), error = function(e) e),
    BPPARAM = BPPARAM
  )
  names(results) <- markers

  # Warn and drop channels where GMM fitting failed
  failed <- vapply(results, inherits, logical(1L), "error")
  if (any(failed)) {
    for (i in which(failed)) {
      warning("Channel '", markers[i], "': bgnorm skipped — ",
              conditionMessage(results[[i]]), call. = FALSE)
    }
    results <- results[!failed]
  }
  if (length(results) == 0L)
    stop("All channels failed bgnorm fitting. Check image quality and cofactor.")

  # Each result is list(adjusted = <vector>, result = <BgnormResult>)
  adj_mat <- vapply(results, function(r) r$adjusted, numeric(nrow(x)))
  colnames(adj_mat) <- names(results)
  br_list <- setNames(lapply(results, function(r) r$result), names(results))
  attr(br_list, "adjusted_matrix") <- adj_mat

  br_list
}
