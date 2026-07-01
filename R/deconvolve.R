#' Background deconvolution for a GMM-fitted intensity vector
#'
#' Implements the conditional-expectation deconvolution from Kharbanda et al.
#' 2025.  Works for both the 3-component (pixel-level) and 2-component
#' (cell-level) models by using the last two components (G-1 and G) as the
#' background and signal pair:
#'
#' \deqn{X_{\text{adj}}(x) = P(C=G \mid x) \Bigl[
#'   (\mu_G - \mu_{G-1}) +
#'   \frac{\sigma_G^2 + \sigma_{G-1}^2 - 2\min(\sigma_{G-1}^2,\sigma_G^2)}
#'        {\sigma_G^2}
#'   (x - \mu_G)
#' \Bigr]}
#'
#' For pixel-level (G=3): component 3 = signal, component 2 = non-specific.
#' For cell-level  (G=2): component 2 = signal, component 1 = non-specific.
#'
#' @param x         Numeric vector of log-transformed intensities.
#' @param means     Ordered GMM means (length G).
#' @param sds       Ordered GMM standard deviations (length G).
#' @param posteriors Matrix (n x G) of posterior probabilities.
#'
#' @return Numeric vector of background-adjusted log-intensities (same length
#'   as \code{x}).
#' @keywords internal
.deconvolve <- function(x, means, sds, posteriors) {
  G      <- length(means)
  mu_bg  <- means[G - 1L];  mu_sig <- means[G]
  v_bg   <- sds[G - 1L]^2;  v_sig  <- sds[G]^2

  var_U  <- v_sig + v_bg - 2 * min(v_bg, v_sig)
  slope  <- var_U / v_sig

  p_sig  <- posteriors[, G]
  pmax(p_sig * ((mu_sig - mu_bg) + slope * (x - mu_sig)), 0)
}

#' Compute the model-based quantile normalisation factor
#'
#' Evaluates the \code{quantile}th quantile of the signal component
#' distribution and maps it through the bgnorm deconvolution formula to
#' obtain a normalisation factor.
#'
#' @param means    Ordered component means (length G; G=3 for pixel, G=2 for
#'   cell).
#' @param sds      Ordered component standard deviations.
#' @param props    Mixing proportions (from the fitted GMM).
#' @param quantile Quantile of the signal component to use (default 0.75).
#'
#' @return Scalar normalisation factor (positive).
#' @keywords internal
.qnorm_factor <- function(means, sds, props, quantile = 0.75) {
  G     <- length(means)
  q_val <- stats::qnorm(quantile, mean = means[G], sd = sds[G])

  # Posterior probability of the signal component at q_val
  w_dens     <- props * stats::dnorm(q_val, means, sds)
  total_dens <- sum(w_dens)
  p_signal   <- if (total_dens > 0) w_dens[G] / total_dens else 0.5

  # Apply the deconvolution formula at the single point q_val
  mu_bg  <- means[G - 1L];  mu_sig <- means[G]
  v_bg   <- sds[G - 1L]^2;  v_sig  <- sds[G]^2
  var_U  <- v_sig + v_bg - 2 * min(v_bg, v_sig)
  q_adj  <- p_signal * ((mu_sig - mu_bg) + (var_U / v_sig) * (q_val - mu_sig))

  if (q_adj <= 0) {
    warning("Quantile normalisation factor <= 0; returning 1.")
    return(1)
  }
  q_adj
}
