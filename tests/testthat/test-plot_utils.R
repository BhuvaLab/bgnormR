# Tests for unified plot_* functions

# -- helpers -----------------------------------------------------------------
.make_spe <- function(n = 60L, mk = c("DAPI", "CD3")) {
    skip_if_not_installed("SpatialExperiment")
    counts <- matrix(
        exp(stats::rnorm(length(mk) * n, log(500), 0.8)),
        nrow = length(mk), ncol = n,
        dimnames = list(mk, paste0("c", seq_len(n)))
    )
    cd <- data.frame(
        x = stats::runif(n), y = stats::runif(n),
        sample_id = rep(c("S1", "S2"), each = n %/% 2L)
    )
    spe <- SpatialExperiment::SpatialExperiment(
        assays     = list(counts = counts),
        colData    = cd,
        sample_id  = "sample_id",
        spatialCoordsNames = c("x", "y")
    )
    bgnorm_sce(spe, assay.type = "counts", name = "bgnorm")
}

# -- plot_qptiff ---------------------------------------------------------------

test_that("plot_qptiff returns ggplot for raw QPTIFFImage", {
    img <- sim_pixel_image()
    p   <- plot_qptiff(img)
    expect_s3_class(p, "ggplot")
})

test_that("plot_qptiff returns ggplot for bgnorm QPTIFFImage", {
    img <- sim_pixel_image()
    res <- bgnorm_pixels(img)
    p   <- plot_qptiff(res)
    expect_s3_class(p, "ggplot")
})

test_that("plot_qptiff returns ggplot for SpatialExperiment", {
    skip_if_not_installed("SpatialExperiment")
    spe <- .make_spe()
    p   <- plot_qptiff(spe, markers = "DAPI")
    expect_s3_class(p, "ggplot")
})

test_that("plot_qptiff errors for non-spatial SE", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:4, 2, 2))
    )
    expect_error(plot_qptiff(se), "SpatialExperiment")
})

# -- plot_pixel_classes --------------------------------------------------------

test_that("plot_pixel_classes returns ggplot for bgnorm QPTIFFImage", {
    img <- sim_pixel_image()
    res <- bgnorm_pixels(img)
    p   <- plot_pixel_classes(res)
    expect_s3_class(p, "ggplot")
})

test_that("plot_pixel_classes errors when QPTIFFImage has no bgnorm results", {
    img <- sim_pixel_image()
    expect_error(plot_pixel_classes(img), "bgnorm results")
})

test_that("plot_pixel_classes returns ggplot for SpatialExperiment", {
    skip_if_not_installed("SpatialExperiment")
    spe <- .make_spe()
    p   <- plot_pixel_classes(spe)
    expect_s3_class(p, "ggplot")
})

test_that("plot_pixel_classes multi-sample SPE includes sample_id in panel labels", {
    skip_if_not_installed("SpatialExperiment")
    spe <- .make_spe(n = 80L)
    p   <- plot_pixel_classes(spe)
    # Build the ggplot to extract the panel data labels
    built <- ggplot2::ggplot_build(p)
    panel_labels <- built$layout$layout$sample_id
    expect_false(is.null(panel_labels))
    expect_true(length(unique(panel_labels)) > 1L)
})

# -- plot_distributions --------------------------------------------------------

test_that("plot_distributions returns ggplot for bgnorm QPTIFFImage", {
    img <- sim_pixel_image()
    res <- bgnorm_pixels(img)
    p   <- plot_distributions(res)
    expect_s3_class(p, "ggplot")
})

test_that("plot_distributions also works with (raw_img, bgnorm_img) two-arg form", {
    img <- sim_pixel_image()
    res <- bgnorm_pixels(img)
    p   <- plot_distributions(img, res)
    expect_s3_class(p, "ggplot")
})

test_that("BgnormResult stores 100-bin histogram of log-intensities", {
    img <- sim_pixel_image()
    res <- bgnorm_pixels(img)
    h   <- bgnorm_results(res)[["Ch1"]]$histogram
    expect_type(h, "list")
    expect_named(h, c("breaks", "density"))
    expect_length(h$density, 100L)
    expect_length(h$breaks, 101L)
    expect_true(all(h$density >= 0))
})

test_that("plot_distributions returns ggplot for SummarizedExperiment", {
    skip_if_not_installed("SpatialExperiment")
    spe <- .make_spe()
    p   <- plot_distributions(spe)
    expect_s3_class(p, "ggplot")
})

test_that("plot_distributions errors when results missing for QPTIFFImage", {
    img <- sim_pixel_image()
    expect_error(plot_distributions(img), "bgnorm results")
})

# -- plot_jsd_heatmap ----------------------------------------------------------

test_that("plot_jsd_heatmap returns ggplot for bgnorm QPTIFFImage", {
    img <- sim_pixel_image()
    res <- bgnorm_pixels(img)
    p   <- plot_jsd_heatmap(res)
    expect_true(inherits(p, "gg"))
})

test_that("plot_jsd_heatmap returns ggplot for SummarizedExperiment", {
    skip_if_not_installed("SpatialExperiment")
    spe <- .make_spe()
    p   <- plot_jsd_heatmap(spe)
    expect_true(inherits(p, "gg"))
})

test_that("plot_jsd_heatmap errors when no bgnorm_results in metadata", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:4, 2, 2))
    )
    expect_error(plot_jsd_heatmap(se), "bgnorm_results")
})

test_that("plot_jsd_heatmap errors when QPTIFFImage has no bgnorm results", {
    img <- sim_pixel_image()
    expect_error(plot_jsd_heatmap(img), "bgnorm results")
})

test_that("plot_jsd_heatmap clusters samples for a list of bgnorm QPTIFFImages", {
    set.seed(1)
    img <- sim_pixel_image(n_markers = 2L, ch_names = c("DAPI", "CD3"))
    res1 <- bgnorm_pixels(img)
    res2 <- bgnorm_pixels(img)
    p    <- plot_jsd_heatmap(list(S1 = res1, S2 = res2))
    expect_true(inherits(p, "gg"))
    # Two samples → y-axis should have 2 levels
    built <- ggplot2::ggplot_build(p)
    expect_equal(length(unique(built$data[[1L]]$y)), 2L)
})
