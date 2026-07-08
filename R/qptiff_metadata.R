## ============================================================
## qptiff_metadata.R
##
## OME-organised metadata model shared by all readers
## (QPTIFF, OME-TIFF, OME-Zarr).  Mirrors the dataclass hierarchy of the
## Python reference (rtubelleza/bioio-tifffile, qptiff_types.py):
##
##   QPTIFFMetadata
##     slide              : slide-constant identity (SlideInfo)
##     images[]           : one scene per distinct image
##       image_info       : optics / acquisition (ImageInfo)
##         camera         : detector settings (CameraInfo)
##         scan_resolution: objective + pixel size (ScanResolutionInfo)
##       channels[]       : per-channel metadata (ChannelInfo)
##       scales[]         : per-pyramid-level metadata (ScaleInfo)
##       raw_xml
##     acquisition_format : brightfield | polaris_scanband | fusion_paged |
##                          ome_tiff | ome_zarr | unknown
##     raw_xml
##
## Consumers use the accessor helpers below rather than hand-writing deep
## list paths, so the schema can evolve without touching call sites.
## ============================================================

#' QPTIFFMetadata: an OME-organised image metadata object
#'
#' An S3-classed nested list returned by \code{\link{metadata}} on a
#' \code{\link{QPTIFFImage}}.  It mirrors the OME hierarchy: a slide-constant
#' \code{slide} block, a list of image \code{images} (each with
#' \code{image_info}, \code{channels}, and \code{scales}), an
#' \code{acquisition_format}, and the raw source \code{raw_xml}.  Prefer the
#' accessor helpers - \code{\link{qpi_format}}, \code{\link{qpi_channels}},
#' \code{\link{qpi_channel_names}}, \code{\link{qpi_pixel_size_um}},
#' \code{\link{qpi_n_levels}}, \code{\link{qpi_primary_scene}},
#' \code{\link{qpi_is_brightfield}}, \code{\link{channel_table}} - over indexing
#' the nested list directly.
#'
#' @return Not applicable; documents the \code{QPTIFFMetadata} structure.
#' @name QPTIFFMetadata
#' @seealso \code{\link{metadata}}, \code{\link{read_qptiff}}
NULL

# ---- Acquisition-format constants (shared) --------------------------------

.QPTIFF_OME_TIFF <- "ome_tiff"
.QPTIFF_OME_ZARR <- "ome_zarr"

# ============================================================
# Constructors
# ============================================================

#' Canonical (empty) SlideInfo block
#' @keywords internal
.new_slide_info <- function(...) {
  base <- list(
    slide_id = NULL, barcode = NULL, study_name = NULL, operator_name = NULL,
    computer_name = NULL, datetime = NULL, acquisition_software = NULL,
    description_version = NULL, instrument_type = NULL, identifier = NULL
  )
  utils::modifyList(base, list(...))
}

#' Canonical (empty) CameraInfo block
#' @keywords internal
.new_camera_info <- function(...) {
  base <- list(camera_name = NULL, camera_type = NULL, gain = NULL,
               bit_depth = NULL, binning = NULL)
  utils::modifyList(base, list(...))
}

#' Canonical (empty) ScanResolutionInfo block
#' @keywords internal
.new_scan_resolution <- function(...) {
  base <- list(magnification = NULL, objective_name = NULL, binning = NULL,
               base_pixel_size_um = NULL)
  utils::modifyList(base, list(...))
}

#' Canonical ImageInfo block (nests camera + scan_resolution)
#' @keywords internal
.new_image_info <- function(..., camera = .new_camera_info(),
                            scan_resolution = .new_scan_resolution()) {
  base <- list(
    image_type = NULL, objective = NULL, bf_lamp_type = NULL,
    scan_profile_name = NULL, scan_mode = NULL, is_tma = NULL,
    opal_kit_type = NULL, xposition_um = NULL, yposition_um = NULL
  )
  out <- utils::modifyList(base, list(...))
  out$camera          <- camera
  out$scan_resolution <- scan_resolution
  out
}

#' Canonical ScaleInfo block (one pyramid level)
#' @keywords internal
.new_scale_info <- function(level, dims = c("y", "x"),
                            downsample_factors = NULL, pixel_sizes_um = NULL) {
  list(level = as.integer(level), dims = dims,
       downsample_factors = downsample_factors, pixel_sizes_um = pixel_sizes_um)
}

#' Canonical scene (one logical image: FullResolution / Label / Macro / ...)
#' @keywords internal
.new_scene <- function(image_info = .new_image_info(), channels = list(),
                       scales = list(), raw_xml = "") {
  list(image_info = image_info, channels = channels, scales = scales,
       raw_xml = raw_xml)
}

#' Canonical QPTIFFMetadata object
#'
#' @param slide              SlideInfo list.
#' @param images             List of scenes (see \code{.new_scene}).
#' @param acquisition_format One of the format constants.
#' @param raw_xml            Raw first-page XML / attrs string.
#' @return A \code{"QPTIFFMetadata"} S3 object.
#' @keywords internal
.new_qptiff_metadata <- function(slide = .new_slide_info(), images = list(),
                                 acquisition_format = .QPTIFF_UNK,
                                 raw_xml = "") {
  structure(
    list(slide = slide, images = images,
         acquisition_format = acquisition_format, raw_xml = raw_xml),
    class = "QPTIFFMetadata"
  )
}

# ============================================================
# Assembly helpers
# ============================================================

#' Build a list of ScaleInfo blocks for an isotropic power-of-two pyramid
#'
#' @param n_levels    Number of pyramid levels (>= 1).
#' @param base_px_um  Full-resolution pixel size in microns, or \code{NULL}.
#' @param base        Downsample base between consecutive levels (default 2).
#' @return List of \code{ScaleInfo} blocks, one per level.
#' @keywords internal
.build_scales <- function(n_levels, base_px_um = NULL, base = 2) {
  n_levels <- max(1L, as.integer(n_levels))
  lapply(seq_len(n_levels), function(l) {
    f  <- base^(l - 1L)
    px <- if (!is.null(base_px_um)) rep(base_px_um * f, 2L) else NULL
    .new_scale_info(level = l, dims = c("y", "x"),
                    downsample_factors = c(f, f), pixel_sizes_um = px)
  })
}

#' Resolve a channel selection to indices + names
#'
#' Shared by all readers: turns a \code{channels} argument (names, 1-based
#' integer indices, or \code{NULL} = all) into aligned \code{ch_idx} /
#' \code{ch_names}, erroring with the available set when a name is not found.
#'
#' @param channels  User channel selection (character / numeric / NULL).
#' @param all_names Character vector of all available channel names.
#' @return List with \code{ch_idx} (integer) and \code{ch_names} (character).
#' @keywords internal
.select_channels <- function(channels, all_names) {
  if (is.null(channels))
    return(list(ch_idx = seq_along(all_names), ch_names = all_names))
  if (is.numeric(channels)) channels <- all_names[as.integer(channels)]
  idx <- match(channels, all_names)
  if (anyNA(idx))
    stop("Channel(s) not found: ", paste(channels[is.na(idx)], collapse = ", "),
         "\nAvailable: ", paste(all_names, collapse = ", "))
  list(ch_idx = as.integer(idx), ch_names = all_names[idx])
}

#' Align a scene's channel metadata to the loaded channels
#'
#' Normalises the primary scene's channel list to \code{length(all_names)}
#' entries (padding with name-only stubs so the list is never shorter than the
#' image), then subsets it to the loaded channels \code{ch_idx}.  Keeps the
#' stored per-channel metadata 1:1 with \code{names(img)}.
#'
#' @param meta      A \code{QPTIFFMetadata} object (single scene).
#' @param all_names Character vector of all channel names (length = n channels).
#' @param ch_idx    Integer indices of the loaded channels.
#' @return \code{meta} with its primary scene's channels aligned and subset.
#' @keywords internal
.subset_scene_channels <- function(meta, all_names, ch_idx) {
  cm <- meta$images[[1L]]$channels
  cm <- lapply(seq_along(all_names), function(i) {
    ch <- if (i <= length(cm)) cm[[i]] else list()
    ch$index          <- ch$index %||% (i - 1L)
    ch$name           <- ch$name  %||% all_names[i]
    ch$is_brightfield <- ch$is_brightfield %||% FALSE
    ch
  })
  meta$images[[1L]]$channels <- cm[ch_idx]
  meta
}

#' Split a flat qpi://vectra vendor map into image- and channel-level fields
#'
#' Vendor keys of the form \code{ch<N>_<field>} go to a per-channel list keyed
#' by \code{<N>}; everything else is an image-level field.  Shared by the
#' OME-TIFF (MapAnnotation) and OME-Zarr (\code{qpi} block) readers.
#'
#' @param kv Named list of string vendor fields.
#' @return List with \code{image} (named list) and \code{channel} (list keyed by
#'   channel-index string).
#' @keywords internal
.split_qpi_map <- function(kv) {
  image <- list(); channel <- list()
  for (k in names(kv)) {
    m <- regmatches(k, regexec("^ch([0-9]+)_(.+)$", k))[[1L]]
    if (length(m) == 3L)
      channel[[m[2L]]] <- c(channel[[m[2L]]], setNames(list(kv[[k]]), m[3L]))
    else
      image[[k]] <- kv[[k]]
  }
  list(image = image, channel = channel)
}

#' Assemble a single-scene QPTIFFMetadata object
#'
#' Convenience wrapper used by every reader: wraps one image scene (the
#' FullResolution scan) into the top-level \code{QPTIFFMetadata} container.
#'
#' @param slide       SlideInfo list.
#' @param image_info  Nested ImageInfo list.
#' @param channels    List of per-channel metadata lists.
#' @param scales      List of ScaleInfo blocks.
#' @param format      Acquisition-format constant.
#' @param raw_xml     Raw description string.
#' @param image_type  Scene image type (default \code{"FullResolution"}).
#' @return A \code{QPTIFFMetadata} object.
#' @keywords internal
.single_scene_metadata <- function(slide, image_info, channels, scales,
                                    format, raw_xml = "",
                                    image_type = "FullResolution") {
  image_info$image_type <- image_info$image_type %||% image_type
  scene <- .new_scene(image_info = image_info, channels = channels,
                      scales = scales, raw_xml = raw_xml)
  .new_qptiff_metadata(slide = slide, images = list(scene),
                       acquisition_format = format, raw_xml = raw_xml)
}

# ============================================================
# Accessors  (mirror the @property helpers in qptiff_types.py)
# ============================================================

#' Acquisition format of a parsed image
#' @param meta A \code{QPTIFFMetadata} object (or legacy metadata list).
#' @return Character scalar, or \code{NULL}.
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' meta <- metadata(read_qptiff(path))
#' qpi_format(meta)
qpi_format <- function(meta) {
  if (is.null(meta)) return(NULL)
  meta$acquisition_format %||% meta$format
}

#' Primary scene of a parsed image (FullResolution if present, else the first)
#' @param meta A \code{QPTIFFMetadata} object.
#' @return A scene list, or \code{NULL} when no scenes are present.
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' scene <- qpi_primary_scene(metadata(read_qptiff(path)))
#' scene$image_info$image_type
qpi_primary_scene <- function(meta) {
  ims <- meta$images
  if (is.null(ims) || length(ims) == 0L) return(NULL)
  types <- vapply(ims, function(im) im$image_info$image_type %||% "", character(1L))
  full  <- which(types == "FullResolution")
  if (length(full) > 0L) ims[[full[1L]]] else ims[[1L]]
}

#' Per-channel metadata of the primary scene
#' @param meta A \code{QPTIFFMetadata} object.
#' @return A list of per-channel metadata lists (possibly empty).
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' length(qpi_channels(metadata(read_qptiff(path))))
qpi_channels <- function(meta) {
  sc <- qpi_primary_scene(meta)
  if (is.null(sc)) return(list())
  sc$channels %||% list()
}

#' Channel names of the primary scene
#' @param meta A \code{QPTIFFMetadata} object.
#' @return Character vector of channel names (possibly empty).
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' qpi_channel_names(metadata(read_qptiff(path)))
qpi_channel_names <- function(meta) {
  ch <- qpi_channels(meta)
  if (length(ch) == 0L) return(character(0L))
  vapply(ch, function(c) c$name %||% NA_character_, character(1L))
}

#' Physical pixel size (microns) of the primary scene, full resolution
#' @param meta A \code{QPTIFFMetadata} object.
#' @return Numeric scalar, or \code{NULL}.
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' qpi_pixel_size_um(metadata(read_qptiff(path)))
qpi_pixel_size_um <- function(meta) {
  sc <- qpi_primary_scene(meta)
  if (is.null(sc)) return(NULL)
  sc$image_info$scan_resolution$base_pixel_size_um
}

#' Number of pyramid levels of the primary scene
#' @param meta A \code{QPTIFFMetadata} object.
#' @return Integer scalar (at least 1).
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' qpi_n_levels(metadata(read_qptiff(path)))
qpi_n_levels <- function(meta) {
  sc <- qpi_primary_scene(meta)
  if (is.null(sc)) return(1L)
  max(1L, length(sc$scales))
}

#' Is the primary scene a brightfield (RGB) image?
#' @param meta A \code{QPTIFFMetadata} object.
#' @return Logical scalar.
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' qpi_is_brightfield(metadata(read_qptiff(path)))
qpi_is_brightfield <- function(meta) {
  fmt <- qpi_format(meta)
  if (!is.null(fmt) && identical(fmt, .QPTIFF_BF)) return(TRUE)
  ch <- qpi_channels(meta)
  any(vapply(ch, function(c) isTRUE(c$is_brightfield), logical(1L)))
}

#' Tidy per-channel metadata table
#'
#' Collapses the primary scene's channel list into a one-row-per-channel
#' \code{data.frame} for quick inspection.  Multi-valued fields such as
#' \code{color_rgb} are rendered as strings.
#'
#' @param meta A \code{QPTIFFMetadata} object.
#' @return A \code{data.frame} with one row per channel; zero rows when the
#'   image has no channel metadata.
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' channel_table(metadata(img))
channel_table <- function(meta) {
  ch <- qpi_channels(meta)
  if (length(ch) == 0L)
    return(data.frame())

  scalar <- function(field, mode = "character") {
    empty <- switch(mode, character = NA_character_, numeric = NA_real_,
                    logical = NA, integer = NA_integer_)
    vapply(ch, function(c) {
      v <- c[[field]]
      if (is.null(v) || length(v) != 1L) return(empty)
      switch(mode, character = as.character(v), numeric = as.numeric(v),
             logical = as.logical(v), integer = as.integer(v))
    }, empty)
  }
  colour <- vapply(ch, function(c) {
    v <- c$color_rgb
    if (is.null(v) || length(v) != 3L) NA_character_
    else paste(v, collapse = ",")
  }, character(1L))

  data.frame(
    index                    = scalar("index", "integer"),
    name                     = scalar("name"),
    fluorophore              = scalar("fluorophore"),
    dye_from_name            = scalar("dye_from_name"),
    exposure_time_us         = scalar("exposure_time_us", "numeric"),
    emission_wavelength_nm   = scalar("emission_wavelength_nm", "numeric"),
    excitation_wavelength_nm = scalar("excitation_wavelength_nm", "numeric"),
    color_rgb                = colour,
    is_brightfield           = scalar("is_brightfield", "logical"),
    stringsAsFactors         = FALSE
  )
}

# ============================================================
# print method
# ============================================================

#' @export
print.QPTIFFMetadata <- function(x, ...) {
  cat("QPTIFFMetadata\n")
  cat("  Format :", qpi_format(x) %||% "unknown", "\n")
  cat("  Scenes :", length(x$images), "\n")
  sl <- x$slide
  if (!is.null(sl$slide_id))  cat("  Slide  :", sl$slide_id, "\n")
  if (!is.null(sl$acquisition_software))
    cat("  Software:", sl$acquisition_software, "\n")
  ch <- qpi_channels(x)
  if (length(ch) > 0L) {
    nms <- qpi_channel_names(x)
    cat("  Channels:", length(ch), "(",
        paste(utils::head(nms, 6L), collapse = ", "),
        if (length(nms) > 6L) "..." else "", ")\n")
  }
  px <- qpi_pixel_size_um(x)
  if (!is.null(px)) cat("  Pixel  :", px, "um\n")
  cat("  Levels :", qpi_n_levels(x), "\n")
  invisible(x)
}
