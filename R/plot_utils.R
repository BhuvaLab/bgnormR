## ============================================================
## plot_utils.R  -  visualisation functions for bgnormR
## ============================================================

# Shared colour palette for GMM components (pixel classes + density lines)
.COMPONENT_COLORS <- c(
  "Background"   = "#2D3436",
  "Non-specific" = "#0984E3",
  "Signal"       = "#D63031"
)
.COMPONENT_LEVELS <- c("Background", "Non-specific", "Signal")

# Qualitative colours for composite multi-channel overlay (up to 10 channels)
.COMPOSITE_COLORS <- c(
  "#00FFFF",  # cyan
  "#FF0000",  # red
  "#00FF00",  # green
  "#FFFF00",  # yellow
  "#FF00FF",  # magenta
  "#FF8000",  # orange
  "#00BFFF",  # deep sky blue
  "#FF69B4",  # hot pink
  "#BF00FF",  # violet
  "#FFFFFF"   # white
)

# ---- internal helpers -------------------------------------------------------

# Downsample a H x W matrix by taking every `factor`-th pixel
.subsample_mat <- function(mat, factor) {
  if (factor <= 1L) return(mat)
  mat[seq(1L, nrow(mat), by = factor),
      seq(1L, ncol(mat), by = factor), drop = FALSE]
}


# ============================================================
# Internal QPTIFFImage (pixel) implementations
# ============================================================

.plot_qptiff_img <- function(img, markers = NULL, resolution = 1L,
                               palette = "magma", scale = "sample", ncol = NULL) {
  scale <- match.arg(scale, c("sample", "marker", "none"))
  MAX_CH      <- 10L
  resolution  <- max(1L, as.integer(resolution))
  d           <- dim(img)
  all_chs     <- dimnames(img)[[3L]] %||% paste0("Ch", seq_len(d[3L]))
  is_adjusted <- !is.null(bgnorm_results(img))

  # Select channels (max MAX_CH; default = first MAX_CH by index)
  if (is.null(markers)) {
    sel_idx <- seq_len(min(d[3L], MAX_CH))
  } else {
    idx  <- match(markers, all_chs)
    miss <- markers[is.na(idx)]
    if (length(miss))
      stop("Markers not found: ", paste(miss, collapse = ", "))
    if (length(idx) > MAX_CH) {
      warning("More than ", MAX_CH, " channels requested; ",
              "displaying first ", MAX_CH, ".")
      idx <- idx[seq_len(MAX_CH)]
    }
    sel_idx <- idx
  }
  n_ch      <- length(sel_idx)
  ch_names  <- all_chs[sel_idx]
  ch_colors <- .COMPOSITE_COLORS[seq_len(n_ch)]

  row_idx <- seq(1L, d[1L], by = resolution)
  col_idx <- seq(1L, d[2L], by = resolution)
  H <- length(row_idx)
  W <- length(col_idx)

  # Pre-compute global bounds once when scaling at the sample level.
  if (scale == "sample") {
    all_vals <- if (.is_lazy_qptiff(img))
      as.vector(as.array(img$.da[row_idx, col_idx, , drop = FALSE]))
    else
      as.vector(unclass(img)[row_idx, col_idx, ])
    if (is_adjusted) all_vals <- 2^all_vals
    lo_global  <- min(all_vals, na.rm = TRUE)
    hi_global  <- stats::quantile(all_vals, 0.999, na.rm = TRUE, names = FALSE)
    rng_global <- hi_global - lo_global
  }

  # Additive RGB composite: each channel contributes intensity x channel_colour.
  rgb_acc <- array(0, dim = c(H, W, 3L))
  for (i in seq_len(n_ch)) {
    mat <- if (.is_lazy_qptiff(img))
      drop(as.array(img$.da[row_idx, col_idx, sel_idx[i], drop = FALSE]))
    else
      unclass(img)[row_idx, col_idx, sel_idx[i]]

    if (is_adjusted) mat <- 2^mat   # invert log2 transform for bgnorm images

    if (scale == "sample") {
      if (rng_global > 0) mat <- pmin((mat - lo_global) / rng_global, 1) else mat[] <- 0
    } else if (scale == "marker") {
      lo  <- min(mat, na.rm = TRUE)
      hi  <- stats::quantile(mat, 0.999, na.rm = TRUE, names = FALSE)
      rng <- hi - lo
      if (rng > 0) mat <- pmin((mat - lo) / rng, 1) else mat[] <- 0
    } else {
      mat <- pmin(pmax(mat, 0), 1)
    }

    col_rgb <- grDevices::col2rgb(ch_colors[i]) / 255  # 3x1 matrix
    rgb_acc[,,1L] <- rgb_acc[,,1L] + mat * col_rgb[1L]
    rgb_acc[,,2L] <- rgb_acc[,,2L] + mat * col_rgb[2L]
    rgb_acc[,,3L] <- rgb_acc[,,3L] + mat * col_rgb[3L]
  }
  rgb_acc <- pmin(rgb_acc, 1)  # clamp: additive overlap saturates toward white

  # Build raster: H x W character matrix of hex colours (row 1 = image top)
  ras <- grDevices::as.raster(
    matrix(
      grDevices::rgb(as.vector(rgb_acc[,,1L]),
                     as.vector(rgb_acc[,,2L]),
                     as.vector(rgb_acc[,,3L])),
      nrow = H, ncol = W
    )
  )

  # Invisible dummy points to anchor the colour legend
  legend_df <- data.frame(
    x      = NA_real_,
    y      = NA_real_,
    marker = factor(ch_names, levels = ch_names)
  )

  ggplot2::ggplot() +
    ggplot2::annotation_raster(ras,
                                xmin = 0.5, xmax = W + 0.5,
                                ymin = 0.5, ymax = H + 0.5,
                                interpolate = FALSE) +
    ggplot2::geom_point(
      data = legend_df,
      ggplot2::aes(x = x, y = y, colour = marker),
      na.rm = TRUE, size = 0
    ) +
    ggplot2::scale_colour_manual(
      values = setNames(ch_colors, ch_names),
      name   = "Marker",
      guide  = ggplot2::guide_legend(
        override.aes = list(size = 4, shape = 15, alpha = 1)
      )
    ) +
    ggplot2::coord_fixed(
      xlim   = c(0.5, W + 0.5),
      ylim   = c(0.5, H + 0.5),
      expand = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = "black", colour = NA),
      panel.background  = ggplot2::element_rect(fill = "black", colour = NA),
      plot.title       = ggplot2::element_text(colour = "white"),
      plot.subtitle       = ggplot2::element_text(colour = "white"),
      plot.caption       = ggplot2::element_text(colour = "white"),
      legend.text       = ggplot2::element_text(colour = "white", size = 9),
      legend.title      = ggplot2::element_text(colour = "white", size = 10,
                                                  face = "bold"),
      legend.background = ggplot2::element_rect(fill = "black", colour = NA),
      legend.key        = ggplot2::element_rect(fill = "black", colour = NA)
    )
}

.plot_pixel_classes_img <- function(results, img, markers = NULL,
                                     resolution = 1L, ncol = NULL) {
  H <- dim(img)[1L]; W <- dim(img)[2L]
  if (!is.null(markers)) results <- results[markers]
  resolution <- max(1L, as.integer(resolution))
  row_idx <- seq(1L, H, by = resolution)
  col_idx <- seq(1L, W, by = resolution)

  df <- do.call(rbind, lapply(names(results), function(nm) {
    res <- results[[nm]]
    adj <- as.vector(as.array(img)[,, nm])
    thr <- res$threshold
    G   <- length(res$parameters$means)
    if (!is.null(thr)) {
      cls <- 1L + (adj > 0) + (adj > thr)
    } else {
      cls <- 1L + (adj > 0)
    }
    labs    <- tail(.COMPONENT_LEVELS, G)
    cls_mat <- .subsample_mat(matrix(cls, nrow = H, ncol = W), resolution)
    data.frame(
      row    = rep(seq_along(row_idx), times = length(col_idx)),
      col    = rep(seq_along(col_idx), each  = length(row_idx)),
      class  = factor(labs[as.vector(cls_mat)], levels = .COMPONENT_LEVELS),
      marker = nm,
      stringsAsFactors = FALSE
    )
  }))
  df$marker <- factor(df$marker, levels = names(results))

  ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = class)) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::facet_wrap(~ marker, ncol = ncol) +
    ggplot2::scale_fill_manual(values   = .COMPONENT_COLORS,
                                name     = "Component",
                                drop     = FALSE,
                                na.value = "black") +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(size = 9, face = "bold",
                                                colour = "white"),
      strip.background = ggplot2::element_rect(fill = "grey20", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "black",  colour = NA),
      panel.background = ggplot2::element_rect(fill = "black",  colour = NA),
      plot.title       = ggplot2::element_text(colour = "white"),
      plot.subtitle       = ggplot2::element_text(colour = "white"),
      plot.caption       = ggplot2::element_text(colour = "white"),
      legend.text      = ggplot2::element_text(colour = "white"),
      legend.title     = ggplot2::element_text(colour = "white"),
      panel.spacing    = grid::unit(3, "pt")
    )
}

.plot_distributions_img <- function(results, markers = NULL, ncol = NULL) {
  avail <- names(results)
  if (!is.null(markers)) avail <- intersect(markers, avail)
  if (length(avail) == 0L)
    stop("No matching channels in bgnorm results.")
  .build_distribution_plot(results[avail], avail, ncol)
}

# ============================================================
# Internal SE / SpatialExperiment implementations
# ============================================================

# Route from SE/SPE for intensity plot - requires SPE for spatial coords
.plot_qptiff_spe <- function(spe, markers = NULL, assay.type = "bgnorm",
                               palette = "magma", point_size = 1,
                               pixels = c(1024L, 1024L), flip_y = TRUE,
                               large_data_threshold = 10000L, ncol = NULL) {
  all_markers <- rownames(spe)
  if (is.null(markers)) markers <- all_markers
  miss <- setdiff(markers, all_markers)
  if (length(miss))
    stop("Markers not found: ", paste(miss, collapse = ", "),
         "\nAvailable: ", paste(head(all_markers, 10L), collapse = ", "))

  base_df <- .spe_base_df(spe)
  mat     <- SummarizedExperiment::assay(spe[markers, , drop = FALSE], assay.type)

  df <- do.call(rbind, lapply(markers, function(mk) {
    v  <- as.numeric(mat[mk, ])
    lo <- min(v, na.rm = TRUE); hi <- max(v, na.rm = TRUE)
    if (hi > lo) v <- (v - lo) / (hi - lo)
    cbind(base_df, value = v, marker = mk, stringsAsFactors = FALSE)
  }))
  df$marker <- factor(df$marker, levels = markers)

  n <- nrow(base_df)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = value)) +
    .spatial_geom(n, point_size, as.integer(pixels), large_data_threshold) +
    .spatial_facet(df, ncol) +
    ggplot2::scale_colour_viridis_c(
      option   = palette, limits = c(0, 1),
      name     = "Intensity\n(scaled)", na.value = "grey30"
    ) +
    ggplot2::coord_equal(expand = FALSE) +
    .spatial_theme()
  if (flip_y) p <- p + ggplot2::scale_y_reverse()
  p
}

# Route from SE/SPE for class plot - requires SPE for spatial coords
.plot_pixel_classes_spe <- function(spe, markers = NULL, point_size = 1,
                                     pixels = c(1024L, 1024L), flip_y = TRUE,
                                     large_data_threshold = 10000L, ncol = NULL) {
  results <- S4Vectors::metadata(spe)$bgnorm_results
  if (is.null(results))
    stop("No bgnorm results in metadata(spe)$bgnorm_results. ",
         "Run bgnorm_sce() first.")

  avail <- names(results)
  if (is.null(markers)) markers <- avail
  miss <- setdiff(markers, avail)
  if (length(miss))
    stop("Markers not in bgnorm results: ", paste(miss, collapse = ", "),
         "\nAvailable: ", paste(head(avail, 10L), collapse = ", "))

  assay_name <- S4Vectors::metadata(spe)$bgnorm_assay %||% "bgnorm"
  adj_mat    <- SummarizedExperiment::assay(spe, assay_name)  # markers x cells

  base_df <- .spe_base_df(spe)

  df <- do.call(rbind, lapply(markers, function(mk) {
    res <- results[[mk]]
    adj <- as.numeric(adj_mat[mk, ])
    thr <- res$threshold
    G   <- length(res$parameters$means)
    if (!is.null(thr)) {
      cls <- 1L + (adj > 0) + (adj > thr)
    } else {
      cls <- 1L + (adj > 0)
    }
    labs <- tail(.COMPONENT_LEVELS, G)
    cbind(base_df,
          class  = factor(labs[cls], levels = .COMPONENT_LEVELS),
          marker = mk,
          stringsAsFactors = FALSE)
  }))
  df$marker <- factor(df$marker, levels = markers)

  n <- nrow(base_df)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = class)) +
    .spatial_geom(n, point_size, as.integer(pixels), large_data_threshold) +
    .spatial_facet(df, ncol) +
    ggplot2::scale_colour_manual(values   = .COMPONENT_COLORS,
                                  name     = "Component",
                                  drop     = FALSE,
                                  na.value = "grey30") +
    ggplot2::coord_equal(expand = FALSE) +
    .spatial_theme()
  if (flip_y) p <- p + ggplot2::scale_y_reverse()
  p
}

# Route from SE/SPE for distributions - uses precomputed histograms from
# metadata(se)$bgnorm_results (stored during bgnorm_sce()).
.plot_distributions_se <- function(se, markers = NULL, ncol = NULL) {
  results <- S4Vectors::metadata(se)$bgnorm_results
  if (is.null(results))
    stop("No bgnorm results in metadata(se)$bgnorm_results. ",
         "Run bgnorm_sce() first.")

  avail <- intersect(rownames(se), names(results))
  if (!is.null(markers)) avail <- intersect(markers, avail)
  if (length(avail) == 0L)
    stop("No matching markers between SE rows and bgnorm results.")
  .build_distribution_plot(results[avail], avail, ncol)
}

# Shared distribution-plot builder.
# results: named list of BgnormResult objects with $histogram (breaks + density).
# markers: character vector of marker names (ordering for facets).
.build_distribution_plot <- function(results, markers, ncol) {
  hist_df <- do.call(rbind, lapply(markers, function(nm) {
    h <- results[[nm]]$histogram
    if (is.null(h))
      stop("BgnormResult for '", nm, "' has no precomputed histogram. ",
           "Re-run bgnorm_pixels() or bgnorm_sce().")
    data.frame(
      xmin    = head(h$breaks, -1L),
      xmax    = tail(h$breaks, -1L),
      density = h$density,
      marker  = nm,
      stringsAsFactors = FALSE
    )
  }))
  hist_df$marker <- factor(hist_df$marker, levels = markers)

  dens_df <- do.call(rbind, lapply(markers, function(nm) {
    p   <- results[[nm]]$parameters
    G   <- length(p$means)
    rng <- range(results[[nm]]$histogram$breaks)
    xs  <- seq(rng[1L], rng[2L], length.out = 512L)
    labs <- tail(.COMPONENT_LEVELS, G)
    do.call(rbind, lapply(seq_len(G), function(k) {
      data.frame(
        x         = xs,
        density   = p$props[k] * stats::dnorm(xs, p$means[k], p$sds[k]),
        component = labs[k],
        marker    = nm,
        stringsAsFactors = FALSE
      )
    }))
  }))
  dens_df$component <- factor(dens_df$component, levels = .COMPONENT_LEVELS)
  dens_df$marker    <- factor(dens_df$marker,    levels = markers)

  ggplot2::ggplot(hist_df) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = density),
      fill = "grey75", colour = NA, alpha = 0.7
    ) +
    ggplot2::geom_line(
      data      = dens_df,
      ggplot2::aes(x = x, y = density, colour = component),
      linewidth = 0.9,
      na.rm     = TRUE,
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~ marker, ncol = ncol, scales = "free") +
    ggplot2::scale_colour_manual(values = .COMPONENT_COLORS,
                                  name   = "Component", drop = FALSE) +
    ggplot2::labs(x = "log2(I/c + 1)", y = "Density") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(size = 9, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom"
    )
}

.plot_jsd_heatmap_impl <- function(results, cluster_rows = TRUE,
                                    cluster_cols = TRUE, show_tissue_positivity = TRUE) {
  is_single <- all(vapply(results, inherits, logical(1L), "BgnormResult"))
  if (is_single) results <- list(sample = results)

  markers    <- Reduce(intersect, lapply(results, names))
  if (length(markers) == 0L)
    stop("No common markers found across samples.")
  sample_nms <- names(results)

  jsd_mat <- do.call(rbind, lapply(results, function(r) {
    vapply(r[markers], `[[`, numeric(1L), "jsd")
  }))
  rownames(jsd_mat) <- sample_nms; colnames(jsd_mat) <- markers

  tp_mat <- do.call(rbind, lapply(results, function(r) {
    vapply(r[markers], function(res) {
      p <- res$parameters$props
      G <- length(p)
      if (G == 3L) p[3L] / (p[2L] + p[3L]) else p[2L]
    }, numeric(1L))
  }))
  rownames(tp_mat) <- sample_nms; colnames(tp_mat) <- markers

  # Replace NA JSD (no-signal channels) with 0 before clustering so that
  # dist() does not produce NA distances and hclust() does not error.
  jsd_clust <- jsd_mat
  jsd_clust[is.na(jsd_clust)] <- 0

  col_ord <- if (cluster_cols && length(markers) > 1L)
    stats::hclust(stats::dist(t(jsd_clust)))$order
  else seq_along(markers)
  row_ord <- if (cluster_rows && length(sample_nms) > 1L)
    stats::hclust(stats::dist(jsd_clust))$order
  else seq_along(sample_nms)

  mk_levels  <- markers[col_ord]
  smp_levels <- sample_nms[row_ord]
  jsd_mat    <- jsd_mat[row_ord, col_ord, drop = FALSE]
  tp_mat     <- tp_mat[row_ord, col_ord, drop = FALSE]

  .mat_to_df <- function(mat, value_col) {
    df <- data.frame(
      sample = rep(rownames(mat), times = ncol(mat)),
      marker = rep(colnames(mat), each  = nrow(mat)),
      value  = as.vector(mat),
      stringsAsFactors = FALSE
    )
    names(df)[3L] <- value_col
    df$sample <- factor(df$sample, levels = smp_levels)
    df$marker <- factor(df$marker, levels = mk_levels)
    df
  }
  jsd_df <- .mat_to_df(jsd_mat, "JSD")

  # Combined overlay: tissue positivity + JSD quality flag for circle colour
  overlay_df <- data.frame(
    sample  = rep(rownames(jsd_mat), times = ncol(jsd_mat)),
    marker  = rep(colnames(jsd_mat), each  = nrow(jsd_mat)),
    tp      = as.vector(tp_mat),
    quality = cut(as.vector(jsd_mat),
                  breaks = c(-Inf, 0.1, 0.2, Inf),
                  labels = c("Low", "Moderate", "Good"),
                  right  = FALSE),
    stringsAsFactors = FALSE
  )
  overlay_df$sample  <- factor(overlay_df$sample,  levels = smp_levels)
  overlay_df$marker  <- factor(overlay_df$marker,   levels = mk_levels)
  overlay_df$quality <- factor(overlay_df$quality,  levels = c("Low", "Moderate", "Good"))

  .hm_theme <- function() {
    list(
      ggplot2::theme_minimal(base_size = 11),
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
        axis.ticks.x     = ggplot2::element_blank(),
        axis.text.y      = ggplot2::element_text(size = 8),
        axis.title       = ggplot2::element_blank(),
        legend.key.width = grid::unit(0.9, "lines")
      )
    )
  }

  p <- ggplot2::ggplot(jsd_df,
                       ggplot2::aes(x = marker, y = sample, fill = JSD)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_viridis_c(option = "F", direction = 1,
                                   limits = c(0, log(2)),
                                   oob    = scales::squish, name = "JSD") +
    ggplot2::scale_x_discrete(position = "bottom") +
    .hm_theme()

  if (show_tissue_positivity) {
    p <- p +
      ggplot2::geom_point(
        data        = overlay_df,
        ggplot2::aes(x = marker, y = sample, size = tp, colour = quality),
        shape       = 19,
        alpha       = 0.85,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_size_area(
        max_size = 8,
        name     = "Tissue\npositivity",
        labels   = function(x) paste0(round(100 * x), "%")
      ) +
      ggplot2::scale_colour_manual(
        values   = c(Low = "red", Moderate = "orange", Good = "white"),
        name     = "JSD quality",
        na.value = "grey60",
        drop     = FALSE
      )
  }

  p
}

# ============================================================
# 1.  plot_qptiff()
# ============================================================

#' Spatial intensity plot for a QPTIFFImage or SpatialExperiment
#'
#' For a \code{\link{QPTIFFImage}}: renders up to 10 channels as an additive
#' colour composite on a black background.  Each channel is assigned a distinct
#' colour; its per-pixel intensity (min-max normalised to [0, 1]) drives the
#' channel's contribution - high-intensity pixels appear fully saturated while
#' low-intensity pixels are transparent (black).  Channels that overlap in
#' space produce mixed additive colours (e.g. cyan + red -> white), matching
#' the standard composite view in FIJI / napari.  When no \code{markers} are
#' provided the first 10 channels (by index) are displayed.
#'
#' For a \code{\link[SpatialExperiment]{SpatialExperiment}}: plots each cell
#' at its spatial coordinates, coloured by its (per-channel scaled) intensity
#' from the requested assay using a viridis colour scale and dark theme.  When
#' the cell count exceeds \code{large_data_threshold} and \pkg{scattermore} is
#' installed, the scatter layer is rasterised automatically.
#'
#' @param x A \code{\link{QPTIFFImage}} \emph{or} a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}}.
#' @param markers Character vector of channel names to display, or \code{NULL}
#'   to use the first 10 channels (for \code{QPTIFFImage}) or all features
#'   (for \code{SpatialExperiment}).  At most 10 channels are shown for
#'   \code{QPTIFFImage} input; excess channels are dropped with a warning.
#' @param resolution Positive integer; pixel downsample factor for
#'   \code{QPTIFFImage} input (ignored for \code{SpatialExperiment}).
#' @param palette Viridis palette name used for \code{SpatialExperiment} input
#'   only.  Default \code{"magma"}.
#' @param assay.type Assay to visualise when \code{x} is a
#'   \code{SpatialExperiment}.  Default \code{"bgnorm"}.
#' @param point_size Base point size for \code{SpatialExperiment} scatter.
#' @param pixels Integer \code{c(width, height)} for \pkg{scattermore}
#'   rasterisation.  Default \code{c(1024L, 1024L)}.
#' @param flip_y Logical; reverse y-axis for \code{SpatialExperiment}?
#'   Default \code{TRUE}.
#' @param large_data_threshold Cell count above which \pkg{scattermore} is used
#'   automatically.  Default \code{10000L}.
#' @param scale Character; intensity scaling for \code{QPTIFFImage} input
#'   (ignored for \code{SpatialExperiment}).  One of:
#'   \describe{
#'     \item{\code{"marker"}}{(default) Per-channel min and 99.9th-percentile; each marker
#'       is stretched to full brightness independently.}
#'     \item{\code{"sample"}}{Global min and 99.9th-percentile computed
#'       across all channels in the image; preserves relative intensities between
#'       markers.}
#'     \item{\code{"none"}}{No scaling; values are clamped to \code{[0, 1]}.}
#'   }
#' @param ncol Ignored for \code{QPTIFFImage} (single composite panel).
#'   Number of facet columns for \code{SpatialExperiment}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot_pixel_classes}}, \code{\link{plot_distributions}},
#'   \code{\link{bgnorm_pixels}}, \code{\link{bgnorm_sce}}
#' @importFrom ggplot2 ggplot aes annotation_raster geom_point geom_raster
#' @importFrom ggplot2 facet_wrap scale_colour_manual scale_fill_viridis_c
#' @importFrom ggplot2 scale_colour_viridis_c guide_legend scale_y_reverse
#' @importFrom ggplot2 coord_equal coord_fixed theme_void theme
#' @importFrom ggplot2 element_text element_blank element_rect labs
#' @importFrom grDevices col2rgb rgb as.raster
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames
#' @importFrom SummarizedExperiment assay colData
#' @importFrom grid unit
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' plot_qptiff(img, markers = c("PanCK", "CD20"))
#'
#' \donttest{
#' # After normalisation: plot background-adjusted intensities
#' res <- bgnorm_pixels(img, sample_prop = 0.1)
#' plot_qptiff(res, markers = c("PanCK", "CD20"), scale = "sample")
#' }
plot_qptiff <- function(x, markers = NULL, resolution = 1L,
                         palette = "magma", assay.type = "bgnorm",
                         scale = c("marker", "sample", "none"),
                         point_size = 1, pixels = c(1024L, 1024L),
                         flip_y = TRUE, large_data_threshold = 10000L,
                         ncol = NULL) {
  scale <- match.arg(scale)
  if (inherits(x, "SpatialExperiment"))
    .plot_qptiff_spe(x, markers = markers, assay.type = assay.type,
                     palette = palette, point_size = point_size,
                     pixels = as.integer(pixels), flip_y = flip_y,
                     large_data_threshold = large_data_threshold, ncol = ncol)
  else if (inherits(x, "SummarizedExperiment"))
    stop("'x' is a SummarizedExperiment but not a SpatialExperiment. ",
         "Spatial coordinates are required for this plot.")
  else if (inherits(x, "QPTIFFImage"))
    .plot_qptiff_img(x, markers = markers, resolution = resolution,
                     palette = palette, scale = scale, ncol = ncol)
  else
    stop("'x' must be a QPTIFFImage or SpatialExperiment.")
}

# ============================================================
# 2.  plot_pixel_classes()
# ============================================================

#' Spatial GMM class assignment plot for a QPTIFFImage or SpatialExperiment
#'
#' Assigns each pixel/cell to its most probable GMM component and renders the
#' result.  For a \code{\link{QPTIFFImage}} produced by
#' \code{\link{bgnorm_pixels}}: pseudocoloured raster (pixel-level data, three
#' components).  For a
#' \code{\link[SpatialExperiment]{SpatialExperiment}}: scatter plot at spatial
#' coordinates using the bgnorm results from
#' \code{metadata(spe)$bgnorm_results} (two components).  Multi-sample
#' \code{SpatialExperiment}s are automatically detected via the
#' \code{sample_id} column in \code{colData} and faceted accordingly.
#'
#' @param x A \code{\link{QPTIFFImage}} returned by \code{\link{bgnorm_pixels}}
#'   (contains both the adjusted pixel intensities and the embedded
#'   \code{\link{BgnormResult}} objects), \emph{or} a
#'   \code{\link[SpatialExperiment]{SpatialExperiment}} with
#'   \code{metadata(x)$bgnorm_results} populated by \code{\link{bgnorm_sce}}.
#' @param markers Character vector of markers to display, or \code{NULL} for
#'   all markers in \code{x}.
#' @param resolution Positive integer; pixel downsample factor for raster
#'   output.  Ignored for \code{SpatialExperiment}.
#' @param point_size Base point size for \code{SpatialExperiment} scatter.
#' @param pixels Integer \code{c(width, height)} for \pkg{scattermore}
#'   rasterisation.  Default \code{c(1024L, 1024L)}.
#' @param flip_y Logical; reverse y-axis for \code{SpatialExperiment}?
#'   Default \code{TRUE}.
#' @param large_data_threshold Cell count threshold for automatic
#'   \pkg{scattermore} use.  Default \code{10000L}.
#' @param ncol Number of facet columns.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot_qptiff}}, \code{\link{plot_distributions}},
#'   \code{\link{bgnorm_pixels}}, \code{\link{bgnorm_sce}}
#' @importFrom ggplot2 ggplot aes geom_raster facet_wrap scale_fill_manual
#' @importFrom ggplot2 scale_colour_manual scale_y_reverse coord_equal
#' @importFrom ggplot2 theme_void theme element_text element_rect
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames
#' @importFrom SummarizedExperiment colData
#' @importFrom grid unit
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' res  <- bgnorm_pixels(img, sample_prop = 0.1)
#' plot_pixel_classes(res, markers = c("PanCK", "CD20"))
plot_pixel_classes <- function(x, markers = NULL,
                                resolution = 1L,
                                point_size = 1,
                                pixels = c(1024L, 1024L),
                                flip_y = TRUE,
                                large_data_threshold = 10000L,
                                ncol = NULL) {
  if (inherits(x, "SpatialExperiment"))
    .plot_pixel_classes_spe(x, markers = markers, point_size = point_size,
                             pixels = as.integer(pixels), flip_y = flip_y,
                             large_data_threshold = large_data_threshold,
                             ncol = ncol)
  else if (inherits(x, "SummarizedExperiment"))
    stop("'x' is a SummarizedExperiment but not a SpatialExperiment. ",
         "Spatial coordinates are required for this plot.")
  else if (inherits(x, "QPTIFFImage")) {
    br <- bgnorm_results(x)
    if (is.null(br))
      stop("'x' is a QPTIFFImage without bgnorm results. ",
           "Run bgnorm_pixels() first and pass its return value.")
    .plot_pixel_classes_img(br, x, markers = markers,
                             resolution = resolution, ncol = ncol)
  } else {
    stop("'x' must be a QPTIFFImage returned by bgnorm_pixels(), ",
         "or a SpatialExperiment.")
  }
}

# ============================================================
# 3.  plot_distributions()
# ============================================================

#' Per-marker intensity distribution with fitted GMM densities
#'
#' Plots a histogram of log2-transformed pixel (or cell) intensities and
#' overlays the component-wise density curves from the fitted GMM.  One facet
#' per marker.  Accepts a \code{\link{QPTIFFImage}} returned by
#' \code{\link{bgnorm_pixels}} (bgnorm results are carried as an attribute) or
#' a \code{\link[SummarizedExperiment]{SummarizedExperiment}} /
#' \code{\link[SpatialExperiment]{SpatialExperiment}} with bgnorm results in
#' \code{metadata(x)$bgnorm_results}.  Multi-sample
#' \code{SpatialExperiment}s are automatically detected via \code{sample_id}
#' in \code{colData} and produce facets by \code{sample_id x marker}.
#'
#' @param x A \code{\link{QPTIFFImage}} returned by \code{\link{bgnorm_pixels}}
#'   \emph{or} a \code{\link[SummarizedExperiment]{SummarizedExperiment}} /
#'   \code{\link[SpatialExperiment]{SpatialExperiment}}.
#' @param results Optionally, the \code{\link{QPTIFFImage}} returned by
#'   \code{\link{bgnorm_pixels}} when \code{x} is the raw (unadjusted)
#'   \code{QPTIFFImage}.  Ignored for \code{SummarizedExperiment} input.
#' @param markers Character vector of markers to display, or \code{NULL} for
#'   all markers.
#' @param ncol   Number of facet columns.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{plot_qptiff}}, \code{\link{plot_pixel_classes}},
#'   \code{\link{bgnorm_pixels}}, \code{\link{bgnorm_sce}}
#' @importFrom ggplot2 ggplot aes geom_rect geom_line facet_wrap
#' @importFrom ggplot2 scale_colour_manual scale_fill_manual labs theme_bw
#' @importFrom ggplot2 theme element_text element_blank element_line
#' @importFrom SummarizedExperiment assay colData
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' res  <- bgnorm_pixels(img, sample_prop = 0.1)
#' plot_distributions(res)
plot_distributions <- function(x, results = NULL, markers = NULL, ncol = NULL) {
  if (inherits(x, "SummarizedExperiment"))
    .plot_distributions_se(x, markers = markers, ncol = ncol)
  else if (inherits(x, "QPTIFFImage")) {
    br <- if (inherits(results, "QPTIFFImage")) {
      bgnorm_results(results) %||%
        stop("'results' QPTIFFImage has no bgnorm results. ",
             "Pass the QPTIFFImage returned by bgnorm_pixels().")
    } else if (is.null(results)) {
      bgnorm_results(x) %||%
        stop("No bgnorm results found. Pass the QPTIFFImage returned by bgnorm_pixels().")
    } else {
      stop("'results' must be the QPTIFFImage returned by bgnorm_pixels().")
    }
    .plot_distributions_img(br, markers = markers, ncol = ncol)
  } else
    stop("'x' must be a QPTIFFImage or SummarizedExperiment.")
}

# ============================================================
# 4.  plot_jsd_heatmap()
# ============================================================

#' QC heatmap of Jensen-Shannon Divergence across markers and samples
#'
#' Displays the per-marker JSD quality metric as a clustered heatmap,
#' optionally annotated with the proportion of the signal GMM component.
#' Accepts a named list of \code{BgnormResult} objects (single sample),
#' a named list of such lists (multiple samples), or a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} /
#' \code{\link[SpatialExperiment]{SpatialExperiment}} with bgnorm results in
#' \code{metadata(results)$bgnorm_results}.
#'
#' @param results A named list of \code{BgnormResult} objects (single sample),
#'   a named list of such lists (multiple samples), a named list of
#'   \code{\link{QPTIFFImage}} objects returned by \code{\link{bgnorm_pixels}}
#'   (one per sample; all must share the same channel names), a single
#'   \code{QPTIFFImage}, or a
#'   \code{SummarizedExperiment} / \code{SpatialExperiment}.
#' @param cluster_rows Logical; cluster samples (rows)?  Default \code{TRUE}.
#' @param cluster_cols Logical; cluster markers (columns)?  Default \code{TRUE}.
#' @param show_tissue_positivity Logical; overlay tissue positivity as circles
#'   on the heatmap?  Circle area is proportional to the tissue positivity
#'   (\eqn{\pi_3 / (\pi_2 + \pi_3)} for three-component models;
#'   \eqn{\pi_2} for two-component models).  Circle colour indicates JSD
#'   quality: red (JSD < 0.1, low), orange (0.1-0.2, moderate), white
#'   (\eqn{\geq} 0.2, good).  Default \code{TRUE}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{bgnorm_pixels}}, \code{\link{bgnorm_sce}}
#' @importFrom ggplot2 ggplot aes geom_tile geom_point scale_fill_viridis_c
#' @importFrom ggplot2 scale_colour_manual scale_size_area scale_x_discrete
#' @importFrom ggplot2 labs theme_minimal theme element_text element_blank
#' @importFrom stats hclust dist
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' res  <- bgnorm_pixels(img, sample_prop = 0.1)
#' plot_jsd_heatmap(res)
#'
#' # Multi-sample comparison (pass a named list)
#' plot_jsd_heatmap(list(sample_A = res, sample_B = res))
plot_jsd_heatmap <- function(results, cluster_rows = TRUE, cluster_cols = TRUE,
                              show_tissue_positivity = TRUE) {
  if (inherits(results, "SummarizedExperiment")) {
    res <- S4Vectors::metadata(results)$bgnorm_results
    if (is.null(res))
      stop("No bgnorm results in metadata(results)$bgnorm_results. ",
           "Run bgnorm_sce() first.")
    results <- res
  } else if (inherits(results, "QPTIFFImage")) {
    res <- bgnorm_results(results)
    if (is.null(res))
      stop("'results' QPTIFFImage has no bgnorm results. ",
           "Pass the QPTIFFImage returned by bgnorm_pixels().")
    results <- res
  } else if (is.list(results) &&
             length(results) > 0L &&
             all(vapply(results, inherits, logical(1L), "QPTIFFImage"))) {
    if (is.null(names(results)))
      stop("The list of QPTIFFImages must be named (one name per sample).")
    br_list <- lapply(results, bgnorm_results)
    nulls   <- vapply(br_list, is.null, logical(1L))
    if (any(nulls))
      stop("All QPTIFFImages must have bgnorm results. ",
           "Run bgnorm_pixels() first on: ",
           paste(names(br_list)[nulls], collapse = ", "))
    results <- br_list
  }
  .plot_jsd_heatmap_impl(results, cluster_rows = cluster_rows,
                          cluster_cols = cluster_cols,
                          show_tissue_positivity = show_tissue_positivity)
}
