# Shared test fixtures for bgnormR unit tests

#' Simulate pixel intensities from the bgnorm generative model
#' Returns raw (non-log) intensities.
sim_pixel_intensities <- function(
    n_bg = 500L, n_ns = 400L, n_sig = 100L,
    mu_bg = 2, mu_ns = 4, mu_sig = 8,
    sd_bg = 0.4, sd_ns = 0.6, sd_sig = 1.0,
    cofactor = 150, seed = 42L) {
  set.seed(seed)
  # Generate log-scale values and transform back
  log_bg  <- rnorm(n_bg,  mu_bg,  sd_bg)
  log_ns  <- rnorm(n_ns,  mu_ns,  sd_ns)
  log_sig <- rnorm(n_sig, mu_sig, sd_sig)
  # Convert from log2(I/150 + 1) scale to raw
  inv_log_transform(c(log_bg, log_ns, log_sig), cofactor = cofactor)
}

#' Simulate cell-level intensities (2-component)
sim_cell_intensities <- function(
    n_ns = 300L, n_sig = 100L,
    mu_ns = 3, mu_sig = 7,
    sd_ns = 0.5, sd_sig = 0.9,
    cofactor = 150, seed = 1L) {
  set.seed(seed)
  log_ns  <- rnorm(n_ns,  mu_ns,  sd_ns)
  log_sig <- rnorm(n_sig, mu_sig, sd_sig)
  inv_log_transform(c(log_ns, log_sig), cofactor = cofactor)
}

#' Wrap pixel intensities into a minimal QPTIFFImage [n x 1 x n_markers]
sim_pixel_image <- function(n_markers = 1L, ch_names = NULL, ...) {
  px <- sim_pixel_intensities(...)
  n  <- length(px)
  if (is.null(ch_names))
    ch_names <- if (n_markers == 1L) "Ch1" else paste0("Ch", seq_len(n_markers))
  mat <- matrix(rep(px, n_markers), nrow = n, ncol = n_markers)
  arr <- array(mat, dim = c(n, 1L, n_markers),
               dimnames = list(NULL, NULL, ch_names))
  as.QPTIFFImage(arr)
}

#' Simulate a multi-marker cell matrix (cells x markers)
sim_cell_matrix <- function(n_cells = 200L, n_markers = 5L, seed = 99L) {
  set.seed(seed)
  mat <- matrix(NA_real_, nrow = n_cells, ncol = n_markers)
  colnames(mat) <- paste0("marker", seq_len(n_markers))
  for (j in seq_len(n_markers)) {
    ns_frac <- sample(c(0.6, 0.7, 0.8), 1L)
    n_ns  <- round(n_cells * ns_frac)
    n_sig <- n_cells - n_ns
    log_vals <- c(rnorm(n_ns, 3, 0.5), rnorm(n_sig, 7, 0.9))
    mat[, j] <- inv_log_transform(log_vals, 150)
  }
  mat
}
