test_that("log_transform is correct", {
  expect_equal(log_transform(0,   150), 0)
  expect_equal(log_transform(150, 150), 1)        # log2(2) = 1
  expect_equal(log_transform(c(0, 150), 150), c(0, 1))
  expect_true(all(log_transform(c(10, 100, 1000)) >= 0))
})

test_that("inv_log_transform inverts log_transform", {
  x <- c(0, 50, 500, 5000, 65535)
  expect_equal(inv_log_transform(log_transform(x)), x, tolerance = 1e-9)
})

test_that("log_transform validates inputs", {
  expect_error(log_transform("a"), "numeric")
  expect_error(log_transform(1:3, cofactor = -1), "cofactor")
  expect_error(log_transform(1:3, cofactor = c(1, 2)), "cofactor")
})

test_that("median_filter_3x3 preserves matrix dimensions", {
  m <- matrix(seq_len(100), 10, 10)
  m2 <- median_filter_3x3(m)
  expect_equal(dim(m), dim(m2))
})

test_that("median_filter_3x3 removes isolated spikes", {
  m <- matrix(0, 7, 7)
  m[4, 4] <- 1000   # isolated spike
  m2 <- median_filter_3x3(m)
  expect_lt(m2[4, 4], 1)  # spike should be largely removed
})

test_that("median_filter_3x3 uses EBImage when available", {
  skip_if_not_installed("EBImage")
  m  <- matrix(0, 9, 9)
  m[5, 5] <- 1000
  m2 <- median_filter_3x3(m)
  expect_equal(dim(m2), dim(m))
  expect_lt(m2[5, 5], 1)
})
