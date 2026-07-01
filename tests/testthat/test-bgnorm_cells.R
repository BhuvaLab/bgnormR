test_that("bgnorm_cells returns BgnormResult with level=cell", {
  x   <- sim_cell_intensities()
  res <- bgnorm_cells(x)
  expect_s3_class(res, "BgnormResult")
  expect_equal(res$level, "cell")
  expect_null(res$adjusted)
  expect_null(res$threshold)   # 2-component: boundary is adj > 0, no stored threshold
})

test_that("bgnorm_cells has 2 components", {
  x   <- sim_cell_intensities()
  res <- bgnorm_cells(x)
  expect_length(res$parameters$means, 2L)
  expect_equal(ncol(res$posteriors), 2L)
})

test_that("bgnorm_cells quantile_norm flag is recorded in result", {
  x    <- sim_cell_intensities()
  r_no <- bgnorm_cells(x, quantile_norm = FALSE)
  r_q  <- bgnorm_cells(x, quantile_norm = TRUE)
  expect_false(r_no$quantile_norm)
  expect_true(r_q$quantile_norm)
})

test_that("bgnorm_sce works on SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  library(SingleCellExperiment)
  counts <- t(sim_cell_matrix(n_cells = 100L, n_markers = 3L))
  sce    <- SingleCellExperiment(assays = list(counts = counts))
  sce2   <- bgnorm_sce(sce, assay.type = "counts", name = "bgnorm")
  expect_true("bgnorm" %in% assayNames(sce2))
  expect_equal(dim(assay(sce2, "bgnorm")), dim(assay(sce2, "counts")))
  expect_true(!is.null(S4Vectors::metadata(sce2)$bgnorm_results))
})

test_that("bgnorm_sce.matrix processes a matrix and returns bgnorm_matrix", {
  mat  <- sim_cell_matrix(n_cells = 80L, n_markers = 4L)
  res  <- bgnorm_sce(mat)
  expect_true(inherits(res, "bgnorm_matrix"))
  expect_equal(dim(res), dim(mat))
})

test_that("bgnorm_cells posteriors have full length when sample_prop < 1", {
  x   <- sim_cell_intensities(n_ns = 200, n_sig = 100)
  res <- bgnorm_cells(x, sample_prop = 0.2)
  expect_equal(nrow(res$posteriors), length(x))
  expect_equal(ncol(res$posteriors), 2L)
})

test_that("bgnorm_sce works on plain SummarizedExperiment", {
  counts <- t(sim_cell_matrix(n_cells = 60L, n_markers = 2L))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  se2 <- bgnorm_sce(se, assay.type = "counts", name = "bgnorm")
  expect_true("bgnorm" %in% SummarizedExperiment::assayNames(se2))
  expect_false(is.null(S4Vectors::metadata(se2)$bgnorm_results))
})

test_that("bgnorm_sce works on SpatialExperiment", {
  skip_if_not_installed("SpatialExperiment")
  counts  <- t(sim_cell_matrix(n_cells = 60L, n_markers = 2L))
  n_cells <- ncol(counts)
  coords  <- data.frame(
    x = runif(n_cells), y = runif(n_cells),
    sample_id = rep("s1", n_cells)
  )
  spe <- SpatialExperiment::SpatialExperiment(
    assays     = list(counts = counts),
    colData    = coords,
    sample_id  = "sample_id",
    spatialCoordsNames = c("x", "y")
  )
  spe2 <- bgnorm_sce(spe, assay.type = "counts", name = "bgnorm")
  expect_true("bgnorm" %in% SummarizedExperiment::assayNames(spe2))
  expect_false(is.null(S4Vectors::metadata(spe2)$bgnorm_results))
})
