#' Log2 transform with cofactor
#'
#' Applies \eqn{I_{log} = \log_2(I / \text{cofactor} + 1)} to stabilise
#' variance in low-intensity fluorescence signals.  The default cofactor of
#' 150 is appropriate for 16-bit images from the Akoya PhenoCycler-Fusion;
#' adjust for other platforms.
#'
#' @param x Numeric vector or matrix of raw intensities.
#' @param cofactor Positive numeric cofactor (default 150 for 16-bit images).
#'
#' @return Numeric vector or matrix of log2-transformed values.
#' @seealso \code{\link{inv_log_transform}}
#' @export
#' @examples
#' x <- c(0, 100, 500, 2000, 10000, 65535)
#' log_transform(x)
log_transform <- function(x, cofactor = 150) {
  stopifnot(is.numeric(x), is.numeric(cofactor), length(cofactor) == 1L,
            cofactor > 0)
  log2(x / cofactor + 1)
}

#' Inverse log2 transform with cofactor
#'
#' Reverses \code{\link{log_transform}}: \eqn{I = (2^{I_{log}} - 1) \times
#' \text{cofactor}}.
#'
#' @param x Numeric vector or matrix of log2-transformed values.
#' @param cofactor Positive numeric cofactor used during \code{log_transform}.
#'
#' @return Numeric vector or matrix on the original intensity scale.
#' @export
#' @examples
#' x <- c(0, 100, 500, 2000, 65535)
#' x_log <- log_transform(x)
#' all.equal(x, inv_log_transform(x_log))
inv_log_transform <- function(x, cofactor = 150) {
  stopifnot(is.numeric(x), is.numeric(cofactor), length(cofactor) == 1L,
            cofactor > 0)
  (2^x - 1) * cofactor
}

#' Apply a 3 x 3 median filter to a matrix
#'
#' Uses \code{EBImage::medianFilter} (C-based, fast) when \pkg{EBImage} is
#' available, otherwise falls back to a pure-R sliding-window implementation.
#'
#' @param img Numeric matrix (height x width).
#'
#' @return Filtered matrix of the same dimensions.
#' @export
#' @examples
#' m <- matrix(runif(100), 10, 10)
#' m_filt <- median_filter_3x3(m)
median_filter_3x3 <- function(img) {
  stopifnot(is.matrix(img), is.numeric(img))

  if (requireNamespace("EBImage", quietly = TRUE)) {
    # EBImage::medianFilter operates on an Image object; w=1 gives a 3x3 kernel
    eb_img <- EBImage::Image(img, colormode = "Grayscale")
    out    <- EBImage::medianFilter(eb_img, size = 1L, cacheSize = 4096L)
    return(matrix(as.numeric(out), nrow = nrow(img), ncol = ncol(img)))
  }

  # Pure-R fallback (border pixels unchanged)
  nr  <- nrow(img); nc <- ncol(img)
  out <- img
  for (r in seq(2L, nr - 1L)) {
    for (col_idx in seq(2L, nc - 1L)) {
      out[r, col_idx] <- stats::median(img[(r - 1L):(r + 1L), (col_idx - 1L):(col_idx + 1L)])
    }
  }
  out
}
