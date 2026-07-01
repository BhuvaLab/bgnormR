#' BgnormResult S3 class
#'
#' @description
#' A lightweight list-based S3 class returned by \code{\link{bgnorm_pixels}}
#' (per channel, via \code{\link{bgnorm_results}}) and \code{\link{bgnorm_cells}}.
#'
#' @field parameters      List with elements \code{means}, \code{sds}, \code{props}.
#' @field posteriors      Matrix (n x G) of posterior probabilities.
#' @field threshold       For the 3-component (pixel-level) model: a single
#'   numeric scalar giving the maximum adjusted intensity of Non-specific pixels.
#'   Combined with the implicit class-1 boundary at 0, this encodes the full
#'   classification: \code{adj == 0} → Background; \code{0 < adj <= threshold}
#'   → Non-specific; \code{adj > threshold} → Signal.
#'   \code{NULL} for 2-component models (cell-level or no-signal fallback),
#'   where the only boundary is 0 (adjusted > 0 → Signal).
#' @field jsd        Jensen-Shannon Divergence QC metric between components 2
#'   and 3.  \code{NA_real_} for no-signal channels.
#' @field level      \code{"pixel"} or \code{"cell"}.
#' @field quantile_norm Logical; whether quantile normalisation was applied.
#' @field histogram  List with \code{$breaks} and \code{$density} vectors from
#'   a pre-computed histogram of log2-transformed intensities.  Used by
#'   \code{\link{plot_distributions}}.
#' @field no_signal  Logical; \code{TRUE} when the 3-component GMM failed and
#'   the channel was fitted with a 2-component fallback, indicating no
#'   detectable biological signal.
#' @field bic Named numeric vector \code{c(G2 = ..., G3 = ...)} with mclust
#'   BIC values (higher is better in mclust's convention) divided by the number
#'   of points used to fit, giving a point-specific estimate comparable across
#'   channels regardless of sample size.  \code{NULL} for cell-level results.
#'
#' @name BgnormResult
NULL

#' @keywords internal
.new_BgnormResult <- function(parameters, posteriors, jsd,
                               level, quantile_norm, no_signal = FALSE,
                               threshold = NULL, histogram = NULL,
                               bic = NULL) {
  structure(
    list(
      parameters    = parameters,
      posteriors    = posteriors,
      threshold     = threshold,
      histogram     = histogram,
      jsd           = jsd,
      bic           = bic,
      level         = level,
      quantile_norm = quantile_norm,
      no_signal     = no_signal
    ),
    class = "BgnormResult"
  )
}

#' @export
print.BgnormResult <- function(x, ...) {
  G    <- length(x$parameters$means)
  props <- x$parameters$props
  cat("BgnormResult (", x$level, "-level)\n", sep = "")
  cat("  n =", nrow(x$posteriors), "\n")
  cat("  Component means:", round(x$parameters$means, 3), "\n")
  cat("  JSD (QC metric):", if (is.na(x$jsd)) "NA" else round(x$jsd, 4), "\n")
  if (!is.null(x$bic)) {
    cat("  BIC (G=2):", round(x$bic["G2"], 2),
        "  BIC (G=3):", round(x$bic["G3"], 2), "\n")
  }
  cat("  No signal detected:", isTRUE(x$no_signal), "\n")
  cat("  Quantile normalised:", x$quantile_norm, "\n")
  if (!isTRUE(x$no_signal)) {
    tp <- if (G == 3L) props[3L] / (props[2L] + props[3L]) else props[2L]
    cat("  Tissue positivity:", paste0(signif(100 * tp, 3), "%"), "\n")
  }
  invisible(x)
}

#' @export
summary.BgnormResult <- function(object, ...) {
  p <- object$parameters
  G <- length(p$means)
  labels <- if (G == 3L) c("Background", "Non-specific", "Signal") else
            c("Non-specific", "Signal")
  df <- data.frame(
    component  = labels,
    mean       = round(p$means, 4),
    sd         = round(p$sds, 4),
    proportion = round(p$props, 4),
    stringsAsFactors = FALSE
  )
  cat("BgnormResult summary (", object$level, "-level)\n", sep = "")
  cat("  JSD:", round(object$jsd, 4), "\n")
  cat("  Quantile normalised:", object$quantile_norm, "\n\n")
  print(df, row.names = FALSE)
  invisible(df)
}
