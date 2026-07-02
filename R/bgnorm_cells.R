#' Cell-level background normalisation (bgnorm)
#'
#' Applies a two-component Gaussian Mixture Model to cell-level aggregated
#' intensities for a single marker.  The two components represent
#' non-specific binding / autofluorescence and biological signal.
#'
#' This is a cell-level approximation of the pixel-level method intended for
#' cases where only cell-summarised intensities (e.g., mean or median
#' intensity per cell) are available.
#'
#' @param x Numeric vector of raw cell-level intensities for one marker.
#'   Must be a plain vector; matrices and arrays are not accepted.
#' @param cofactor Positive numeric cofactor for log2 transform (default 150).
#' @param quantile_norm Logical; apply model-based quantile normalisation?
#' @param quantile Quantile of the signal component for normalisation
#'   (default 0.75).
#' @param sample_prop Numeric in (0, 1]; proportion of non-zero observations
#'   used to fit the GMM.  Default \code{1} (use all).
#' @param ... Additional arguments forwarded to the internal GMM fitter.
#'
#' @return A \code{BgnormResult} object.
#' @seealso \code{\link{bgnorm_pixels}}, \code{\link{bgnorm_sce}}
#' @export
#' @examples
#' set.seed(3)
#' # Simulate cell intensities: nonspecific + signal
#' x <- c(exp(rnorm(300, log(200), 0.5)), exp(rnorm(100, log(1500), 0.8)))
#' res <- bgnorm_cells(x)
#' print(res)
bgnorm_cells <- function(x, cofactor = 150, quantile_norm = FALSE,
                          quantile = 0.75, sample_prop = 1, ...) {
  .bgnorm_cell_impl(x, cofactor = cofactor, quantile_norm = quantile_norm,
                    quantile = quantile, sample_prop = sample_prop, ...)$result
}

# Private implementation that returns list(adjusted, result) for use by
# .bgnorm_columns() so the adjusted matrix can be stored in the SCE assay.
.bgnorm_cell_impl <- function(x, cofactor = 150, quantile_norm = FALSE,
                               quantile = 0.75, sample_prop = 1, ...) {
  if (!is.numeric(x) || !is.null(dim(x)))
    stop("'x' must be a numeric vector.")
  if (length(x) < 20L)
    stop("Too few observations for bgnorm (need >= 20).")
  stopifnot(is.numeric(cofactor), cofactor > 0)

  # 1. Log2-transform
  x_log <- log_transform(x, cofactor = cofactor)

  # 2. Fit 2-component GMM (zeros filtered and optionally subsampled inside)
  gmm <- .fit_gmm(x_log, n_components = 2L, sample_prop = sample_prop, ...)

  # 3. Compute posteriors for the FULL x_log using the fitted parameters.
  #    (gmm$posteriors covers only the fitting subset; we need all cells.)
  full_post <- .gmm_posteriors(x_log, gmm$means, gmm$props, gmm$sds)

  # 4. Order components: lower mean = nonspecific, higher = signal
  ord <- order(gmm$means)
  params <- list(
    means      = gmm$means[ord],
    sds        = gmm$sds[ord],
    props      = gmm$props[ord],
    posteriors = full_post[, ord, drop = FALSE]
  )

  # 5. Deconvolve
  x_adj <- .deconvolve(x_log, params$means, params$sds, params$posteriors)

  # 6. Quantile normalisation
  if (quantile_norm) {
    q_factor <- .qnorm_factor(params$means, params$sds, params$props,
                               quantile = quantile)
    x_adj <- x_adj / q_factor
  }

  # JSD between component 1 (nonspecific) and component 2 (signal)
  jsd <- jsd_gaussians(params$means[1L], params$sds[1L],
                        params$means[2L], params$sds[2L])

  list(
    adjusted = x_adj,
    result   = .new_BgnormResult(
      parameters    = list(means = params$means, sds = params$sds,
                           props = params$props),
      n             = length(x_log),
      threshold     = NULL,   # 2-component: class boundary is adj > 0
      histogram     = .compute_histogram(x_log),
      jsd           = jsd,
      level         = "cell",
      quantile_norm = quantile_norm
    )
  )
}

#' Apply cell-level bgnorm to a SingleCellExperiment or matrix
#'
#' Runs \code{\link{bgnorm_cells}} independently for each marker (column) and
#' returns an updated object with adjusted intensities in a new assay.
#'
#' @param x         A \code{SummarizedExperiment} (including
#'   \code{SingleCellExperiment} and \code{SpatialExperiment} subclasses) or a
#'   numeric matrix (cells \eqn{\times} markers).
#' @param assay.type Character; name of the assay to normalise. Default
#'   \code{"counts"}.
#' @param name      Character; name for the output assay (SCE) or returned
#'   matrix attribute.  Default \code{"bgnorm"}.
#' @param cofactor  Cofactor for log2 transform (default 150).
#' @param quantile_norm Logical; apply bgnormQ?
#' @param quantile  Quantile for normalisation.
#' @param BPPARAM   A \code{\link[BiocParallel]{BiocParallelParam}} instance.
#' @param ...       Additional arguments forwarded to \code{bgnorm_cells}.
#'
#' @return The input object with a new assay (\code{name}) holding the
#'   adjusted intensities.  The per-marker \code{BgnormResult} list is stored
#'   in \code{metadata(x)$bgnorm_results}.
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom S4Vectors metadata<-
#' @importFrom BiocParallel bplapply SerialParam
#' @export
#' @examples
#' library(SummarizedExperiment)
#' set.seed(5)
#' counts <- matrix(
#'   exp(rnorm(2000, log(300), 0.9)),
#'   nrow = 10, ncol = 200,
#'   dimnames = list(paste0("marker", seq_len(10)), paste0("cell", seq_len(200)))
#' )
#' se <- SummarizedExperiment(assays = list(counts = counts))
#' se <- bgnorm_sce(se, assay.type = "counts", name = "bgnorm")
#' assayNames(se)
bgnorm_sce <- function(x, assay.type = "counts", name = "bgnorm",
                        cofactor = 150, quantile_norm = FALSE,
                        quantile = 0.75,
                        BPPARAM = BiocParallel::SerialParam(), ...) {
  UseMethod("bgnorm_sce")
}

#' @export
bgnorm_sce.SummarizedExperiment <- function(x, assay.type = "counts",
                                             name = "bgnorm",
                                             cofactor = 150,
                                             quantile_norm = FALSE,
                                             quantile = 0.75,
                                             BPPARAM = BiocParallel::SerialParam(),
                                             ...) {
  mat <- t(SummarizedExperiment::assay(x, assay.type))  # cells x markers
  res <- bgnorm_sce.matrix(mat,
                            name          = name,
                            cofactor      = cofactor,
                            quantile_norm = quantile_norm,
                            quantile      = quantile,
                            BPPARAM       = BPPARAM, ...)
  adj <- attr(res, "adjusted_matrix")              # cells x markers
  SummarizedExperiment::assay(x, name) <- t(adj)  # markers x cells
  S4Vectors::metadata(x)$bgnorm_results <- attr(res, "results_list")
  S4Vectors::metadata(x)$bgnorm_assay   <- name
  x
}

#' @export
bgnorm_sce.matrix <- function(x, assay.type = NULL, name = "bgnorm",
                               cofactor = 150, quantile_norm = FALSE,
                               quantile = 0.75,
                               BPPARAM = BiocParallel::SerialParam(), ...) {
  if (!is.null(assay.type))
    warning("'assay.type' is ignored when 'x' is a matrix; ",
            "supply the assay matrix directly.", call. = FALSE)
  results <- .bgnorm_columns(
    x       = x,
    fn      = function(col) .bgnorm_cell_impl(col, cofactor = cofactor,
                                               quantile_norm = quantile_norm,
                                               quantile = quantile, ...),
    BPPARAM = BPPARAM
  )
  adj_mat <- attr(results, "adjusted_matrix")
  attr(results, "adjusted_matrix") <- NULL

  structure(adj_mat,
            results_list    = results,
            adjusted_matrix = adj_mat,
            class = c("bgnorm_matrix", "matrix", "array"))
}
