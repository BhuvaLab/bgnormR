test_that("jsd_gaussians is zero for identical distributions", {
  jsd <- jsd_gaussians(mu1 = 1, sd1 = 1, mu2 = 1, sd2 = 1)
  expect_equal(jsd, 0, tolerance = 1e-4)
})

test_that("jsd_gaussians increases with separation", {
  j1 <- jsd_gaussians(0, 1, 1, 1)   # small separation
  j2 <- jsd_gaussians(0, 1, 5, 1)   # large separation
  expect_gt(j2, j1)
})

test_that("jsd_gaussians is bounded by 1 (log2 base)", {
  j <- jsd_gaussians(0, 0.1, 100, 0.1)  # essentially disjoint
  expect_lte(j, 1 + 1e-4)
})

test_that("jsd_gaussians validates inputs", {
  expect_error(jsd_gaussians(0, -1, 1, 1), "sd1")
  expect_error(jsd_gaussians(0,  1, 1, 0), "sd2")
})

test_that("jsd_qc extracts metric from a BgnormResult", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  j   <- jsd_qc(bgnorm_results(res)[["Ch1"]])
  expect_length(j, 1L)
  expect_gte(j, 0)
})

test_that("qc_summary accepts a QPTIFFImage returned by bgnorm_pixels", {
  img <- sim_pixel_image(n_markers = 2L, ch_names = c("CD3", "CD8"),
                          n_bg = 300L, n_ns = 200L, n_sig = 100L)
  res <- bgnorm_pixels(img)
  df  <- qc_summary(res)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("name", "jsd", "prop_signal"))
  expect_equal(nrow(df), 2L)
  expect_setequal(df$name, c("CD3", "CD8"))
})

test_that("qc_summary accepts a named list of BgnormResult", {
  img1 <- sim_pixel_image(ch_names = "CD3", seed = 11L)
  img2 <- sim_pixel_image(ch_names = "CD8", seed = 22L)
  r1   <- bgnorm_pixels(img1)
  r2   <- bgnorm_pixels(img2)
  df   <- qc_summary(list(CD3 = bgnorm_results(r1)[["CD3"]],
                           CD8 = bgnorm_results(r2)[["CD8"]]))
  expect_s3_class(df, "data.frame")
  expect_named(df, c("name", "jsd", "prop_signal"))
  expect_equal(nrow(df), 2L)
})
