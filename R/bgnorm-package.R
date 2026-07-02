#' bgnormR: Background Normalisation for Multiplex Spatial Proteomics
#'
#' Implements \emph{bgnorm}, a generative statistical framework for
#' background correction, normalisation, and quality control in multiplex
#' spatial proteomics (Kharbanda et al. 2025).
#'
#' The framework models log-transformed fluorescence intensities as a
#' three-component Gaussian mixture representing:
#' \enumerate{
#'   \item Background (empty space / instrument noise)
#'   \item Non-specific binding and autofluorescence
#'   \item True biological signal
#' }
#'
#' Core functions:
#' \itemize{
#'   \item \code{\link{bgnorm_pixels}} - pixel-level background correction
#'   \item \code{\link{bgnorm_cells}} - cell-level background correction
#'   \item \code{\link{bgnorm_sce}} - normalise a SummarizedExperiment / SpatialExperiment
#'   \item \code{\link{jsd_qc}} - Jensen-Shannon Divergence quality metric
#'   \item \code{\link{read_qptiff}} - Akoya PhenoCycler-Fusion image reader
#'   \item \code{\link{plot_qptiff}}, \code{\link{plot_pixel_classes}},
#'         \code{\link{plot_distributions}}, \code{\link{plot_jsd_heatmap}} - plotting utilities
#' }
#'
#' @docType package
#' @name bgnormR-package
#' @aliases bgnormR
#' @importFrom stats setNames
#' @importFrom utils tail
#' @importFrom graphics hist
"_PACKAGE"

# Column names referenced via ggplot2 NSE (aes()) rather than as R symbols.
utils::globalVariables(c(
  "JSD", "component", "density", "marker", "quality", "tp", "value",
  "x", "xmax", "xmin", "y"
))
