## ============================================================
## plot_spatial.R  -  internal helpers for SPE scatter plots
##
## These helpers are called by plot_qptiff(), plot_pixel_classes(), and
## plot_distributions() when their first argument is a SpatialExperiment.
## ============================================================

# Internal: return the appropriate point layer depending on dataset size.
.spatial_geom <- function(n, point_size, pixels, threshold) {
  if (n > threshold && requireNamespace("scattermore", quietly = TRUE)) {
    scattermore::geom_scattermore(pointsize = point_size, pixels = pixels)
  } else {
    ggplot2::geom_point(size = point_size / 3, stroke = 0)
  }
}

# Internal: dark theme for spatial class plots
.spatial_theme <- function() {
  list(
    ggplot2::theme_void(),
    ggplot2::theme(
      strip.text       = ggplot2::element_text(size = 9, face = "bold",
                                               colour = "white"),
      strip.background = ggplot2::element_rect(fill = "grey20", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "black",  colour = NA),
      panel.background = ggplot2::element_rect(fill = "black",  colour = NA),
      legend.text      = ggplot2::element_text(colour = "white"),
      legend.title     = ggplot2::element_text(colour = "white"),
      panel.spacing    = grid::unit(3, "pt")
    )
  )
}

# Internal: extract spatial coordinates and sample_id
.spe_base_df <- function(spe) {
  coords    <- SpatialExperiment::spatialCoords(spe)
  coord_nms <- SpatialExperiment::spatialCoordsNames(spe)
  sid       <- as.character(SummarizedExperiment::colData(spe)[["sample_id"]])
  data.frame(
    x = coords[, coord_nms[1L]], y = coords[, coord_nms[2L]],
    sample_id = sid, stringsAsFactors = FALSE
  )
}

# Internal: facet by marker, adding sample_id when multiple samples exist
.spatial_facet <- function(df, ncol) {
  if (length(unique(df$sample_id)) > 1L)
    ggplot2::facet_wrap(~ sample_id + marker, ncol = ncol)
  else
    ggplot2::facet_wrap(~ marker, ncol = ncol)
}

# ============================================================
# plot_spatial_classes()  - GMM class scatter plot for SpatialExperiment
# ============================================================
