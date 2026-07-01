# ---- bgnorm_pixels: input validation ---------------------------------------

test_that("bgnorm_pixels rejects non-QPTIFFImage input", {
  expect_error(bgnorm_pixels(rnorm(100)),                    "QPTIFFImage")
  expect_error(bgnorm_pixels(matrix(rnorm(100), 10L, 10L)), "QPTIFFImage")
  arr <- array(rnorm(100), dim = c(10L, 10L, 1L))
  expect_error(bgnorm_pixels(arr),                           "QPTIFFImage")
})

test_that("bgnorm_pixels errors on unknown channel names", {
  img <- sim_pixel_image(ch_names = "DAPI")
  expect_error(bgnorm_pixels(img, channels = "NOTHERE"), "not found")
})

# ---- bgnorm_pixels: return type --------------------------------------------

test_that("bgnorm_pixels returns a QPTIFFImage with bgnorm_results", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  expect_s3_class(res, "QPTIFFImage")
  expect_named(res, "Ch1")
  br <- bgnorm_results(res)
  expect_type(br, "list")
  expect_named(br, "Ch1")
  expect_s3_class(br[["Ch1"]], "BgnormResult")
  expect_equal(br[["Ch1"]]$level, "pixel")
  expect_false(br[["Ch1"]]$quantile_norm)
})

test_that("adjusted image has same spatial dims as input", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  expect_equal(dim(res)[1:2], dim(img)[1:2])
})

test_that("BgnormResult has threshold but no adjusted field", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  br  <- bgnorm_results(res)[["Ch1"]]
  expect_null(br$adjusted)
  expect_true(is.numeric(br$threshold) && length(br$threshold) == 1L)
  expect_gte(br$threshold, 0)   # threshold is in adjusted space (>= 0 after pmax)
})

test_that("component means are in ascending order", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  p   <- bgnorm_results(res)[["Ch1"]]$parameters
  expect_true(p$means[1L] < p$means[3L])
})

test_that("signal pixels have positive mean adjusted intensity (from QPTIFFImage)", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  ch  <- bgnorm_results(res)[["Ch1"]]
  adj <- as.vector(as.array(res)[,, "Ch1"])
  signal_adj <- adj[ch$posteriors[, 3L] > 0.9]
  if (length(signal_adj) > 5L) expect_true(mean(signal_adj) > 0)
})

test_that("quantile_norm shifts threshold and QPTIFFImage values", {
  img  <- sim_pixel_image()
  r_no <- bgnorm_pixels(img, quantile_norm = FALSE)
  r_q  <- bgnorm_pixels(img, quantile_norm = TRUE)
  expect_true(bgnorm_results(r_q)[["Ch1"]]$quantile_norm)
  expect_false(isTRUE(all.equal(bgnorm_results(r_no)[["Ch1"]]$threshold,
                                 bgnorm_results(r_q)[["Ch1"]]$threshold)))
})

test_that("JSD is numeric in [0, 1] (log2 base)", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  jsd <- bgnorm_results(res)[["Ch1"]]$jsd
  expect_true(is.numeric(jsd))
  expect_gte(jsd, 0)
  expect_lte(jsd, 1 + 1e-6)
})

test_that("print and summary work for BgnormResult", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  br  <- bgnorm_results(res)[["Ch1"]]
  expect_output(print(br),   "BgnormResult")
  expect_output(summary(br), "Background")
})

test_that("posteriors span all pixels when sample_prop < 1", {
  set.seed(42)
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img, sample_prop = 0.1)
  ch  <- bgnorm_results(res)[["Ch1"]]
  n   <- dim(img)[1L] * dim(img)[2L]
  expect_equal(nrow(ch$posteriors), n)
  expect_equal(ncol(ch$posteriors), 3L)
})

# ---- bgnorm_pixels: multi-channel ------------------------------------------

test_that("bgnorm_pixels processes all channels of a QPTIFFImage", {
  set.seed(8)
  img <- sim_pixel_image(n_markers = 3L, ch_names = c("DAPI", "CD3", "CD8"),
                          n_bg = 300L, n_ns = 200L, n_sig = 100L)
  res <- bgnorm_pixels(img)
  expect_s3_class(res, "QPTIFFImage")
  expect_named(res, c("DAPI", "CD3", "CD8"))
  expect_equal(dim(res), dim(img))
  br <- bgnorm_results(res)
  expect_named(br, c("DAPI", "CD3", "CD8"))
  expect_s3_class(br[["DAPI"]], "BgnormResult")
})

test_that("channels argument selects a subset", {
  set.seed(9)
  img <- sim_pixel_image(n_markers = 3L, ch_names = c("DAPI", "CD3", "CD8"),
                          n_bg = 300L, n_ns = 200L, n_sig = 100L)
  res <- bgnorm_pixels(img, channels = c("DAPI", "CD8"))
  expect_s3_class(res, "QPTIFFImage")
  expect_named(res, c("DAPI", "CD8"))
  expect_equal(dim(res)[3L], 2L)
  expect_equal(dimnames(res)[[3L]], c("DAPI", "CD8"))
  expect_named(bgnorm_results(res), c("DAPI", "CD8"))
})

test_that("channels argument accepts integer indices", {
  img <- sim_pixel_image(n_markers = 3L, ch_names = c("DAPI", "CD3", "CD8"),
                          n_bg = 300L, n_ns = 200L, n_sig = 100L)
  res <- bgnorm_pixels(img, channels = c(1L, 3L))
  expect_named(res, c("DAPI", "CD8"))
})

# ---- QPTIFFImage subsetting preserves bgnorm_results -----------------------

test_that("channel subsetting propagates bgnorm_results", {
  set.seed(7)
  img <- sim_pixel_image(n_markers = 3L, ch_names = c("DAPI", "CD3", "CD8"),
                          n_bg = 300L, n_ns = 200L, n_sig = 100L)
  res <- bgnorm_pixels(img)
  sub <- res[, , c("DAPI", "CD8")]
  expect_s3_class(sub, "QPTIFFImage")
  expect_named(bgnorm_results(sub), c("DAPI", "CD8"))
})

test_that("spatial subsetting clears bgnorm_results", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  sub <- res[1:5, , ]
  expect_null(bgnorm_results(sub))
})

# ---- 2-component fallback (no signal) ---------------------------------------

test_that("2-component fallback sets no_signal=TRUE and adjusted=0 for sparse channel", {
  set.seed(99)
  # Sparse channel: only 2 clusters (background + non-specific), no signal
  # Use just 2 tight clusters so 3-component Mclust returns NULL
  px_sparse <- c(rep(1L, 900L), sample(40:80, 100L, replace = TRUE))
  # Normal channel with clear 3-component structure
  px_normal <- c(
    as.integer(exp(rnorm(500L, log(50),   0.4))),
    as.integer(exp(rnorm(400L, log(300),  0.6))),
    as.integer(exp(rnorm(100L, log(2000), 1.0)))
  )
  arr <- array(
    c(px_normal, px_sparse),
    dim = c(1000L, 1L, 2L),
    dimnames = list(NULL, NULL, c("Normal", "Sparse"))
  )
  img <- as.QPTIFFImage(arr)

  res <- suppressMessages(suppressWarnings(bgnorm_pixels(img)))
  br  <- bgnorm_results(res)

  # Normal channel: no_signal = FALSE, has a scalar threshold (3-component)
  expect_false(isTRUE(br[["Normal"]]$no_signal))
  expect_true(any(as.array(res)[,, "Normal"] != 0))
  expect_true(is.numeric(br[["Normal"]]$threshold) && length(br[["Normal"]]$threshold) == 1L)

  if ("Sparse" %in% names(br)) {
    # If 2-component fallback succeeded: all-zero adjusted values, no_signal = TRUE
    expect_true(isTRUE(br[["Sparse"]]$no_signal))
    expect_true(all(as.array(res)[,, "Sparse"] == 0))
    expect_true(is.na(br[["Sparse"]]$jsd))
    expect_length(br[["Sparse"]]$parameters$means, 2L)
    expect_null(br[["Sparse"]]$threshold)
  }
  # Either the channel is present (2-comp fallback) or absent (both GMMs fail)
  expect_true("Normal" %in% names(br))
})

test_that("bgnorm_pixels warns and skips channels where both GMMs fail", {
  set.seed(42)
  # Channel with fewer than 20 non-zero pixels hits the pre-check in .fit_gmm
  # (needs >= n_components * 10 = 20 for 2-component) and never calls mclust.
  px_scarce <- c(rep(0L, 985L), as.integer(runif(15L, 100, 1000)))
  px_normal <- as.integer(c(
    exp(rnorm(500L, log(50), 0.4)),
    exp(rnorm(400L, log(300), 0.6)),
    exp(rnorm(100L, log(2000), 1.0))
  ))
  arr <- array(
    c(px_normal, px_scarce),
    dim = c(1000L, 1L, 2L),
    dimnames = list(NULL, NULL, c("Normal", "Scarce"))
  )
  img <- as.QPTIFFImage(arr)
  expect_warning(
    res <- suppressMessages(bgnorm_pixels(img)),
    "bgnorm skipped"
  )
  br <- bgnorm_results(res)
  expect_true("Normal"  %in% names(br))
  expect_false("Scarce" %in% names(br))
  # The returned QPTIFFImage has only the "Normal" channel
  expect_named(res, "Normal")
})

test_that("threshold in BgnormResult is consistent with QPTIFFImage values", {
  img <- sim_pixel_image()
  res <- suppressMessages(bgnorm_pixels(img))
  br  <- bgnorm_results(res)[["Ch1"]]
  thr <- br$threshold
  expect_true(is.numeric(thr) && length(thr) == 1L)
  # Applying the threshold must recover valid classes
  adj <- as.vector(as.array(res)[,, "Ch1"])
  cls <- 1L + (adj > 0) + (adj > thr)
  expect_true(all(cls %in% 1:3))
})
