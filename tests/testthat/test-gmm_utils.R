test_that(".fit_gmm returns correct structure", {
  set.seed(1)
  x   <- c(rnorm(200, 2, 0.3), rnorm(150, 4, 0.5), rnorm(50, 7, 0.8))
  gmm <- bgnormR:::.fit_gmm(x, n_components = 3L)
  expect_named(gmm, c("means", "sds", "props", "posteriors"))
  expect_length(gmm$means, 3L)
  # nrow(posteriors) == number of non-zero finite observations used
  expect_true(nrow(gmm$posteriors) <= length(x))
  expect_equal(ncol(gmm$posteriors), 3L)
  expect_equal(sum(gmm$props), 1, tolerance = 1e-6)
})

test_that(".fit_gmm filters out zeros before fitting", {
  set.seed(2)
  x_clean   <- c(rnorm(200, 2, 0.3), rnorm(150, 4, 0.5), rnorm(50, 7, 0.8))
  x_with_zeros <- c(x_clean, rep(0, 100))
  gmm <- bgnormR:::.fit_gmm(x_with_zeros, n_components = 3L)
  # Posteriors only for non-zero rows
  expect_equal(nrow(gmm$posteriors), length(x_clean))
})

test_that(".fit_gmm respects sample_prop when dataset exceeds threshold", {
  set.seed(3)
  x <- c(rnorm(300, 2, 0.3), rnorm(200, 5, 0.5), rnorm(100, 9, 0.8))
  gmm_full    <- bgnormR:::.fit_gmm(x, n_components = 3L, sample_prop = 1)
  # Force sampling by setting threshold below n_nonzero (600)
  gmm_sampled <- bgnormR:::.fit_gmm(x, n_components = 3L, sample_prop = 0.2,
                                     .sample_threshold = 300)
  expect_lt(nrow(gmm_sampled$posteriors), nrow(gmm_full$posteriors))
  expect_length(gmm_sampled$means, 3L)
  expect_equal(sum(gmm_sampled$props), 1, tolerance = 1e-6)
})

test_that(".fit_gmm emits message and uses all data when below threshold", {
  set.seed(3)
  x <- c(rnorm(300, 2, 0.3), rnorm(200, 5, 0.5), rnorm(100, 9, 0.8))
  expect_message(
    gmm <- bgnormR:::.fit_gmm(x, n_components = 3L, sample_prop = 0.1),
    "Sampling not required"
  )
  # All 600 non-zero observations should be used
  expect_equal(nrow(gmm$posteriors), 600L)
})

test_that(".fit_gmm fails gracefully on tiny inputs", {
  expect_error(bgnormR:::.fit_gmm(rnorm(5), 3L), "Too few")
})

test_that(".resolve_components orders by mean", {
  means <- c(5, 2, 8)
  sds   <- c(0.5, 0.4, 1.0)
  props <- c(0.3, 0.5, 0.2)
  post  <- matrix(rep(props, each = 10), nrow = 10)
  res   <- bgnormR:::.resolve_components(means, sds, props, post)
  expect_equal(res$means[1L], 2)
  expect_equal(res$means[3L], 8)
})

test_that(".resolve_components swaps when right tail rule applies", {
  # Component with lower mean but much larger SD has heavier right tail
  means <- c(1, 3, 4)      # ord = 1,2,3
  sds   <- c(0.2, 2.0, 0.3)  # component 2 (mean=3) has larger SD
  # q75 for (3, 2.0) = 3 + 0.674*2 = 4.35
  # q75 for (4, 0.3) = 4 + 0.674*0.3 = 4.20
  # So component with mean=3 has heavier right tail → swap → signal=mean3 component
  props <- c(0.5, 0.3, 0.2)
  post  <- matrix(rep(props, each = 10), nrow = 10)
  res   <- bgnormR:::.resolve_components(means, sds, props, post)
  # After swap: component 3 should be the one with mean=3 (heavier tail)
  expect_equal(res$means[3L], 3)
  expect_equal(res$means[2L], 4)
})

test_that(".classify_pixels applies threshold corrections for 3-component model", {
  means <- c(1, 3, 6)   # background=1, non-specific=3, signal=6
  # Construct posteriors where argmax favours component 3 (signal) for all pixels
  post <- matrix(c(0.05, 0.05, 0.90), nrow = 3, ncol = 3, byrow = TRUE)

  # x[1] = 0.5 < μ₁ = 1 → should become 1 (background) via pass 1 then pass 2
  # x[2] = 2.0 between μ₁ and μ₂ → should become 2 (non-specific) via pass 1
  # x[3] = 5.0 > μ₂ → argmax stays 3 (signal)
  x   <- c(0.5, 2.0, 5.0)
  cls <- bgnormR:::.classify_pixels(x, means, post)
  expect_equal(cls, c(1L, 2L, 3L))
})

test_that(".classify_pixels demotes cls==2 below μ₁ to background", {
  means <- c(2, 5, 9)
  # posteriors where argmax is component 2 for first pixel, component 3 otherwise
  post <- matrix(
    c(0.1, 0.8, 0.1,
      0.1, 0.1, 0.8),
    nrow = 2, byrow = TRUE
  )
  x   <- c(1.0, 7.0)   # first pixel below μ₁; second above μ₂
  cls <- bgnormR:::.classify_pixels(x, means, post)
  expect_equal(cls[1L], 1L)
  expect_equal(cls[2L], 3L)
})

test_that(".classify_pixels does not apply threshold correction for 2-component model", {
  means <- c(2, 7)
  post  <- matrix(c(0.1, 0.9), nrow = 1)
  cls   <- bgnormR:::.classify_pixels(0.5, means, post)
  # No correction for 2-component; argmax wins
  expect_equal(cls, 2L)
})

test_that(".compute_threshold returns max adjusted value of class-2 pixels", {
  x_adj <- c(0.0, 0.0, 0.5, 1.2, 1.8)
  cls   <- c(1L,  1L,  2L,  2L,  3L)
  thr   <- bgnormR:::.compute_threshold(x_adj, cls)
  expect_equal(thr, 1.2)
})

test_that(".compute_threshold returns 0 when no class-2 pixels", {
  x_adj <- c(0.0, 0.0, 1.5, 2.0)
  cls   <- c(1L,  1L,  3L,  3L)   # no class-2 pixels
  thr   <- bgnormR:::.compute_threshold(x_adj, cls)
  expect_equal(thr, 0)
})
