## ============================================================
## read_ome_zarr.R
##
## OME-NGFF (OME-Zarr) reader.  Loads a multiscale OME-Zarr store written in
## the layout produced by the Python reference (bioio-tifffile qptiff_zarr.py):
## root attributes carry `ome.multiscales` (dataset paths + scale transforms +
## axes) and `ome.omero.channels` (labels + colours), plus an optional top-level
## `qpi` vendor block.  Pixel data is read through Rarr; large stores are
## exposed lazily via a DelayedArray backed by OMEZarrArraySeed.
##
## Requires the Rarr package (Suggests).
## ============================================================

#' @importFrom DelayedArray extract_array type
NULL

# ---- Detection -----------------------------------------------------------

#' Does a path point at an OME-Zarr / zarr store?
#' @keywords internal
.is_ome_zarr_store <- function(path) {
  if (grepl("\\.zarr/?$", path, ignore.case = TRUE)) return(TRUE)
  if (dir.exists(path)) {
    fs <- list.files(path, all.files = TRUE)
    if (any(c("zarr.json", ".zattrs", ".zgroup") %in% fs)) return(TRUE)
  }
  FALSE
}

# ============================================================
# OMEZarrArraySeed  (lazy on-disk zarr access)
# ============================================================

#' Seed class for lazy on-disk OME-Zarr access
#'
#' Satisfies the \pkg{DelayedArray} seed contract (\code{dim}, \code{dimnames},
#' \code{type}, \code{extract_array}) for a single pyramid level of an OME-Zarr
#' store, reading requested blocks through \code{Rarr::read_zarr_array}.
#' Zarr axes (typically \code{c, y, x}) are permuted to the QPTIFFImage
#' convention \code{[y, x, c]} = \code{[H, W, C]}.  Do not construct directly;
#' use \code{\link{read_qptiff}} on a \code{.ome.zarr} store.
#'
#' @slot array_path    Absolute path to the level array inside the store.
#' @slot .dim          Integer \code{c(H, W, C)} in output order.
#' @slot .dimnames     List of dimnames; element 3 holds channel names.
#' @slot dtype         \code{"integer"} or \code{"double"}.
#' @slot axis_names    Character vector of zarr axis names (store order).
#' @slot zarr_channels Integer; 1-based zarr channel index per output channel.
#' @slot divisor       Numeric; value pixels are divided by in \code{"double"} mode.
#' @slot metadata      A \code{QPTIFFMetadata} object.
#'
#' @exportClass OMEZarrArraySeed
setClass("OMEZarrArraySeed", representation(
  array_path    = "character",
  .dim          = "integer",
  .dimnames     = "list",
  dtype         = "character",
  axis_names    = "character",
  zarr_channels = "integer",
  divisor       = "numeric",
  metadata      = "ANY"
))

#' @rdname OMEZarrArraySeed-class
#' @export
setMethod("dim", "OMEZarrArraySeed", function(x) x@.dim)

#' @rdname OMEZarrArraySeed-class
#' @export
setMethod("dimnames", "OMEZarrArraySeed", function(x) {
  nms <- x@.dimnames
  if (all(vapply(nms, is.null, logical(1L)))) NULL else nms
})

#' @rdname OMEZarrArraySeed-class
#' @export
setMethod("type", "OMEZarrArraySeed", function(x) x@dtype)

#' Seed-contract methods for OMEZarrArraySeed
#'
#' @param x     An \code{OMEZarrArraySeed}.
#' @param index List of length 3; each element is \code{NULL} (all indices) or
#'   an integer vector of 1-based indices for that dimension.
#' @return \code{extract_array} returns an ordinary array; \code{dim} an integer
#'   vector; \code{dimnames} a list; \code{type} a character scalar.
#' @rdname OMEZarrArraySeed-class
#' @export
setMethod("extract_array", "OMEZarrArraySeed", function(x, index) {
  an   <- x@axis_names
  H    <- x@.dim[1L]; W <- x@.dim[2L]; Cout <- x@.dim[3L]
  i_rows <- index[[1L]] %||% seq_len(H)
  i_cols <- index[[2L]] %||% seq_len(W)
  i_chn  <- index[[3L]] %||% seq_len(Cout)

  empty_val <- if (x@dtype == "integer") NA_integer_ else NA_real_
  if (length(i_rows) == 0L || length(i_cols) == 0L || length(i_chn) == 0L)
    return(array(empty_val, c(length(i_rows), length(i_cols), length(i_chn))))

  zc <- x@zarr_channels[i_chn]

  # Build the per-axis index in the store's own axis order.
  zidx <- vector("list", length(an))
  for (a in seq_along(an))
    zidx[[a]] <- switch(an[a], y = i_rows, x = i_cols, c = zc, 1L)

  raw <- Rarr::read_zarr_array(x@array_path, index = zidx)

  # Permute (store order) -> (y, x, c); index has collapsed z/t to length 1.
  perm <- match(c("y", "x", "c"), an)
  rest <- setdiff(seq_along(an), perm)
  raw  <- aperm(raw, c(perm, rest))
  d    <- dim(raw)
  dim(raw) <- d[seq_len(3L)]

  if (x@dtype == "integer") {
    storage.mode(raw) <- "integer"
  } else {
    raw <- as.numeric(raw) / x@divisor
    dim(raw) <- d[seq_len(3L)]
  }
  raw
})

# ============================================================
# Attribute parsing
# ============================================================

# NGFF root attrs: v0.5 nests everything under "ome"; v0.4 puts multiscales /
# omero at the root.  Return the block that actually holds `multiscales`.
.ngff_ome_block <- function(attrs) {
  if (!is.null(attrs$ome) && !is.null(attrs$ome$multiscales)) return(attrs$ome)
  attrs
}

#' Build a QPTIFFMetadata object from OME-Zarr root attributes
#' @keywords internal
.zarr_to_metadata <- function(attrs, channel_names, colours, pixel_size_um,
                              n_levels) {
  qpi         <- attrs$qpi %||% list()
  split       <- .split_qpi_map(qpi)   # image-level vs per-channel (ch<N>_*)
  qpi_image   <- split$image
  qpi_channel <- split$channel

  slide <- .new_slide_info(
    slide_id             = qpi_image$slide_id,
    barcode              = qpi_image$barcode,
    study_name           = qpi_image$study_name,
    computer_name        = qpi_image$computer_name,
    datetime             = qpi_image$datetime,
    acquisition_software = qpi_image$acquisition_software,
    description_version  = qpi_image$description_version,
    identifier           = qpi_image$identifier
  )
  image_info <- .new_image_info(
    image_type      = qpi_image$image_type,
    objective       = qpi_image$objective,
    scan_mode       = qpi_image$scan_mode,
    scan_resolution = .new_scan_resolution(base_pixel_size_um = pixel_size_um)
  )

  channels <- lapply(seq_along(channel_names), function(i) {
    out <- list(index = i - 1L, name = channel_names[i], is_brightfield = FALSE)
    if (!is.null(colours) && i <= length(colours) && !is.null(colours[[i]]))
      out$color_rgb <- colours[[i]]
    qc <- qpi_channel[[as.character(i - 1L)]]
    if (!is.null(qc)) out <- utils::modifyList(out, .coerce_qpi_channel(qc))
    out
  })

  .single_scene_metadata(
    slide      = slide,
    image_info = image_info,
    channels   = channels,
    scales     = .build_scales(n_levels, pixel_size_um),
    format     = .QPTIFF_OME_ZARR
  )
}

# Convert a 6-char hex colour ("RRGGBB") to c(r, g, b).
.hex_to_rgb <- function(hex) {
  if (is.null(hex) || is.na(hex)) return(NULL)
  hex <- sub("^#", "", hex)
  if (nchar(hex) < 6L) return(NULL)
  as.integer(strtoi(substring(hex, c(1L, 3L, 5L), c(2L, 4L, 6L)), 16L))
}

# ============================================================
# Top-level OME-Zarr reader
# ============================================================

#' Read an OME-Zarr (OME-NGFF) store into a QPTIFFImage
#' @inheritParams read_qptiff
#' @keywords internal
.read_ome_zarr <- function(path, channels = NULL, level = 1L,
                           as_integer = TRUE, lazy = FALSE) {
  if (!requireNamespace("Rarr", quietly = TRUE))
    stop("Reading OME-Zarr requires the 'Rarr' package.\n",
         "  Install with: BiocManager::install('Rarr')")

  message("Reading OME-Zarr metadata ...")
  attrs <- tryCatch(Rarr::read_zarr_attributes(path),
                    error = function(e)
                      stop("Could not read zarr attributes from: ", path,
                           "\n  ", conditionMessage(e)))
  ome <- .ngff_ome_block(attrs)

  ms <- ome$multiscales
  if (is.null(ms) || length(ms) == 0L)
    stop("No 'multiscales' metadata found; not a recognised OME-Zarr store: ", path)
  ms <- ms[[1L]]

  # --- Axes (store order, e.g. c/y/x) -----------------------------------
  axis_names <- tolower(vapply(ms$axes, function(a) a$name %||% "", character(1L)))

  # --- Datasets / level selection ---------------------------------------
  datasets <- ms$datasets
  n_levels <- length(datasets)
  if (level > n_levels)
    stop("Requested level ", level, " exceeds available pyramid levels (",
         n_levels, ").")
  ds <- datasets[[level]]
  array_path <- file.path(path, ds$path)

  # Physical pixel size from the level-0 scale transform (y, x axes).
  scale0 <- tryCatch(datasets[[1L]]$coordinateTransformations[[1L]]$scale,
                     error = function(e) NULL)
  yx <- match(c("y", "x"), axis_names)
  pixel_size_um <- if (!is.null(scale0) && !is.na(yx[1L])) scale0[[yx[1L]]] else NULL

  # --- Array shape / dtype ----------------------------------------------
  ov <- Rarr::zarr_overview(array_path, as_data_frame = TRUE)
  shape <- as.integer(ov$dim[[1L]])
  dtype <- as.character(ov$data_type[[1L]])

  cpos <- match("c", axis_names)
  ypos <- match("y", axis_names)
  xpos <- match("x", axis_names)
  if (is.na(ypos) || is.na(xpos))
    stop("OME-Zarr store is missing spatial (y/x) axes.")
  n_ch <- if (!is.na(cpos)) shape[cpos] else 1L
  H    <- shape[ypos]; W <- shape[xpos]

  # --- Channel names / colours ------------------------------------------
  omero    <- ome$omero %||% attrs$omero
  ch_meta  <- if (!is.null(omero)) omero$channels else NULL
  all_names <- if (!is.null(ch_meta))
    vapply(ch_meta, function(c) c$label %||% NA_character_, character(1L))
  else rep(NA_character_, n_ch)
  if (length(all_names) < n_ch)
    all_names <- c(all_names, rep(NA_character_, n_ch - length(all_names)))
  all_names <- all_names[seq_len(n_ch)]
  miss <- is.na(all_names) | !nzchar(all_names)
  all_names[miss] <- paste0("Channel_", which(miss) - 1L)

  colours <- if (!is.null(ch_meta))
    lapply(ch_meta, function(c) .hex_to_rgb(c$color)) else NULL

  # --- Channel selection -------------------------------------------------
  sel      <- .select_channels(channels, all_names)
  ch_idx   <- sel$ch_idx
  ch_names <- sel$ch_names

  meta <- .zarr_to_metadata(attrs, all_names, colours, pixel_size_um, n_levels)
  # Keep the stored channel metadata aligned 1:1 with the loaded channels.
  meta <- .subset_scene_channels(meta, all_names, ch_idx)

  # --- Build seed / array ------------------------------------------------
  # 8-bit stores normalise by 255; everything else by 65535 (double mode only).
  divisor <- if (grepl("8$", dtype)) 255 else 65535
  seed <- new("OMEZarrArraySeed",
    array_path    = normalizePath(array_path),
    .dim          = as.integer(c(H, W, length(ch_idx))),
    .dimnames     = list(NULL, NULL, ch_names),
    dtype         = if (as_integer) "integer" else "double",
    axis_names    = axis_names,
    zarr_channels = as.integer(ch_idx),
    divisor       = divisor,
    metadata      = meta
  )

  da <- DelayedArray::DelayedArray(seed)
  if (lazy) return(.new_QPTIFFImage_lazy(da, metadata = meta))

  message("Loading ", length(ch_idx), " channel(s) ...")
  arr <- as.array(da)
  dimnames(arr) <- list(NULL, NULL, ch_names)
  .new_QPTIFFImage(arr, meta)
}
