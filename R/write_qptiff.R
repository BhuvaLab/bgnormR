#' Write a QPTIFFImage to a TIFF file
#'
#' Writes a \code{\link{QPTIFFImage}} as a multi-page 16-bit grayscale
#' OME-TIFF (one page per channel).  The first page carries an OME-XML
#' \code{ImageDescription} (tag 270) that declares the pages as \emph{channels}
#' (\code{SizeC = nChannels}, \code{SizeZ = SizeT = 1}), so viewers such as
#' QuPath / Bio-Formats interpret each page as a channel rather than a
#' timepoint.  Per-channel bgnorm results are embedded so the file is
#' self-documenting and round-trips through \code{\link{read_qptiff}}.
#'
#' @section Intensity transform:
#' \describe{
#'   \item{bgnorm-adjusted images}{The QPTIFFImage stores background-adjusted
#'     log\eqn{_2}-intensities.  These are inverted with \eqn{2^x} before
#'     writing so the output values are in a linear intensity scale.}
#'   \item{Raw images (no bgnorm results)}{Written as-is; values are assumed
#'     to be in \code{[0, 65535]}.}
#' }
#' All output values are rounded to the nearest integer and clamped to the
#' \code{[0, 65535]} range before being written as 16-bit unsigned integers.
#'
#' The native, Java-free TIFF/QPTIFF reader and writer implemented here was
#' translated from the
#' \href{https://github.com/rtubelleza/bioio-tifffile/tree/feature/read-qptiffs-rich}{bioio-tifffile fork}
#' by Rafael Tubelleza.
#'
#' @section Metadata:
#' The first page's \code{ImageDescription} tag (TIFF tag 270) is an OME-XML
#' document.  Its \code{<Pixels>} element declares \code{SizeC = nChannels}
#' with one \code{<Channel Name="...">} per channel and an explicit
#' \code{<TiffData>} page-to-channel mapping, which is what makes downstream
#' viewers read the pages as channels.  Channel 0's bgnorm metadata is appended
#' to the OME root as no-namespace elements (ignored by OME readers but read
#' back by \code{\link{read_qptiff}}); pages 1..n each carry a minimal
#' \code{PerkinElmerQPI} block:
#' \preformatted{
#' <PerkinElmerQPI>
#'   <Biomarker>CD20</Biomarker>
#'   <transform>2^x</transform>
#'   <bgnorm>{"level":"pixel","jsd":0.35,...}</bgnorm>
#' </PerkinElmerQPI>
#' }
#' The \code{<Biomarker>} element is the channel name; \code{<bgnorm>} holds a
#' JSON object with GMM parameters, JSD, threshold, and normalisation flags
#' (only present for bgnorm-adjusted images).  If \pkg{jsonlite} is not
#' installed, the \code{<bgnorm>} element is omitted.
#'
#' @param x             A \code{\link{QPTIFFImage}} (eager or lazy).
#' @param path          Output file path (character scalar).  The directory
#'   must exist.
#' @return \code{path}, invisibly.
#'
#' @seealso \code{\link{read_qptiff}}, \code{\link{bgnorm_pixels}}
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' res  <- bgnorm_pixels(img, sample_prop = 0.1)
#' out  <- file.path(tempdir(), "PA_HNC_bgnorm.tif")
#' write_qptiff(res, out)
#'
#' # Round-trip: channel names and dimensions are preserved
#' img2 <- read_qptiff(out)
#' names(img2)
#' dim(img2)
write_qptiff <- function(x, path) {
  if (!inherits(x, "QPTIFFImage"))
    stop("'x' must be a QPTIFFImage.")
  if (!is.character(path) || length(path) != 1L || !nzchar(path))
    stop("'path' must be a single non-empty character string.")

  chs  <- names(x)
  d    <- dim(x)
  H    <- d[1L]; W <- d[2L]; C <- d[3L]
  br   <- bgnorm_results(x)
  is_adjusted <- !is.null(br)

  transform <- if (is_adjusted) "2^x" else "none"

  # Rich per-channel + slide/image metadata for the OME-XML header, aligned
  # positionally to the loaded channels by name.
  meta    <- attr(x, "metadata")
  ch_meta <- .match_channel_meta(meta, chs)

  # Extract 3D array (materialises lazy images)
  arr <- as.array.QPTIFFImage(x)

  # Build uint16 matrices and per-channel description lists
  mats       <- vector("list", C)
  desc_lists <- vector("list", C)

  for (k in seq_len(C)) {
    ch  <- chs[k]
    mat <- arr[, , k]

    if (is_adjusted) mat <- 2^mat

    # Clamp and round to [0, 65535]
    mat <- pmin(pmax(round(mat), 0), 65535)
    mats[[k]] <- mat

    # Build description
    desc_list <- list(channel = ch, transform = transform)
    if (is_adjusted && !is.null(br[[ch]])) {
      r <- br[[ch]]
      desc_list$bgnorm <- list(
        level         = r$level,
        jsd           = r$jsd,
        no_signal     = isTRUE(r$no_signal),
        quantile_norm = isTRUE(r$quantile_norm),
        threshold     = r$threshold,
        parameters    = list(
          means = r$parameters$means,
          sds   = r$parameters$sds,
          props = r$parameters$props
        )
      )
    }
    desc_lists[[k]] <- desc_list
  }

  # Pages 1..n carry a PerkinElmerQPI block; the first page instead carries an
  # OME-XML description declaring the pages as channels (so QuPath / Bio-Formats
  # read them as channels, not timepoints).  Channel 0's bgnorm metadata is
  # embedded in the OME root as no-namespace elements so read_qptiff still
  # recovers it.
  descs <- vapply(desc_lists, .desc_to_string, character(1L))
  descs[1L] <- .build_ome_description(chs, H, W, C, desc_lists[[1L]],
                                      meta = meta, ch_meta = ch_meta)

  .write_tiff_binary(path, mats, descs, H, W)
  invisible(path)
}

# ---- Internal: build the OME-XML description for the first page -----------
#
# Declares one Image with SizeC = number of channels (SizeZ = SizeT = 1) and an
# explicit page->channel (TiffData) mapping, so OME-aware readers interpret each
# TIFF page as a channel.  Channel 0's Biomarker / transform / bgnorm elements
# are appended to the OME root with an empty namespace: OME readers ignore
# foreign-namespace elements, while read_qptiff's direct-child XPath lookups
# still find them, preserving the round-trip for channel 0.

.build_ome_description <- function(chs, H, W, C, ch0_desc_list,
                                   meta = NULL, ch_meta = NULL) {
  ns <- "http://www.openmicroscopy.org/Schemas/OME/2016-06"
  if (is.null(ch_meta)) ch_meta <- vector("list", C)

  px_um <- qpi_pixel_size_um(meta)
  phys  <- if (!is.null(px_um))
    sprintf(' PhysicalSizeX="%s" PhysicalSizeXUnit="\u00b5m" PhysicalSizeY="%s" PhysicalSizeYUnit="\u00b5m"',
            px_um, px_um)
  else ""

  # Per-channel <Channel> elements: Name + any rich fields we hold.
  channels <- paste0(vapply(seq_len(C), function(k)
    .ome_channel_element(k - 1L, chs[k], ch_meta[[k]]),
    character(1L)), collapse = "")

  # Per-channel <Plane> elements carry exposure time (Z=T=0).
  planes <- paste0(vapply(seq_len(C), function(k)
    .ome_plane_element(k - 1L, ch_meta[[k]]),
    character(1L)), collapse = "")

  tiffdata <- paste0(vapply(seq_len(C), function(k)
    sprintf('<TiffData FirstC="%d" FirstZ="0" FirstT="0" IFD="%d" PlaneCount="1"/>',
            k - 1L, k - 1L),
    character(1L)), collapse = "")

  # PerkinElmer-specific fields with no OME equivalent -> qpi://vectra map.
  annotation <- .qpi_map_annotation(meta, ch_meta, C)
  ann_ref    <- if (nzchar(annotation)) '<AnnotationRef ID="Annotation:0"/>' else ""

  img_name <- .xml_escape(meta$images[[1L]]$image_info$scan_profile_name %||% "bgnormR")

  paste0(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<OME xmlns="', ns, '" ',
    'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
    'xsi:schemaLocation="', ns, ' ', ns, '/ome.xsd">',
    '<Image ID="Image:0" Name="', img_name, '">',
    '<Pixels ID="Pixels:0" DimensionOrder="XYCZT" Type="uint16" ',
    'SizeX="', W, '" SizeY="', H, '" SizeC="', C, '" SizeZ="1" SizeT="1" ',
    'Interleaved="false" BigEndian="false" SignificantBits="16"', phys, '>',
    channels, planes, tiffdata,
    '</Pixels>', ann_ref, '</Image>',
    annotation,
    .desc_body(ch0_desc_list, ns_reset = TRUE),
    '</OME>'
  )
}

# ---- Internal: OME <Channel> element from per-channel metadata ------------

.ome_channel_element <- function(index, name, cm) {
  cm <- cm %||% list()
  attrs <- sprintf(' ID="Channel:0:%d" Name="%s"', index, .xml_escape(name))
  if (!is.null(cm$fluorophore))
    attrs <- paste0(attrs, sprintf(' Fluor="%s"', .xml_escape(cm$fluorophore)))
  col <- .rgb_to_ome_color(cm$color_rgb)
  if (!is.null(col))
    attrs <- paste0(attrs, sprintf(' Color="%s"', col))
  if (!is.null(cm$emission_wavelength_nm))
    attrs <- paste0(attrs, sprintf(' EmissionWavelength="%s" EmissionWavelengthUnit="nm"',
                                   cm$emission_wavelength_nm))
  if (!is.null(cm$excitation_wavelength_nm))
    attrs <- paste0(attrs, sprintf(' ExcitationWavelength="%s" ExcitationWavelengthUnit="nm"',
                                   cm$excitation_wavelength_nm))
  paste0('<Channel', attrs, ' SamplesPerPixel="1"/>')
}

# ---- Internal: OME <Plane> element (exposure) -----------------------------

.ome_plane_element <- function(index, cm) {
  cm <- cm %||% list()
  if (is.null(cm$exposure_time_us)) return("")
  sprintf('<Plane TheZ="0" TheT="0" TheC="%d" ExposureTime="%s" ExposureTimeUnit="\u00b5s"/>',
          index, cm$exposure_time_us)
}

# ---- Internal: encode c(r, g, b) as an OME signed-int32 RGBA colour --------

.rgb_to_ome_color <- function(rgb) {
  if (is.null(rgb) || length(rgb) != 3L || anyNA(rgb)) return(NULL)
  u <- rgb[1L] * 2^24 + rgb[2L] * 2^16 + rgb[3L] * 2^8 + 255   # alpha = 255
  if (u >= 2^31) u <- u - 2^32                                  # to signed int32
  format(as.integer(u), scientific = FALSE)
}

# ---- Internal: match loaded channel names to metadata channels ------------

.match_channel_meta <- function(meta, chs) {
  if (is.null(meta)) return(vector("list", length(chs)))
  cm <- qpi_channels(meta)
  if (length(cm) == 0L) return(vector("list", length(chs)))
  cm_names <- vapply(cm, function(c) c$name %||% "", character(1L))
  # Consume each metadata entry at most once so repeated channel names (e.g.
  # a marker imaged in several cycles) map to distinct metadata entries in
  # order rather than all collapsing onto the first match.
  used <- logical(length(cm))
  out  <- vector("list", length(chs))
  for (j in seq_along(chs)) {
    i <- which(cm_names == chs[j] & !used)
    if (length(i) == 0L) next          # out[[j]] stays NULL -> no metadata
    used[i[1L]] <- TRUE
    out[[j]]    <- cm[[i[1L]]]
  }
  out
}

# ---- Internal: qpi://vectra MapAnnotation (round-trips PerkinElmer fields) --
#
# Mirrors qptiff_ome.py: image-level keys plus per-channel ch<N>_* keys for
# fields that have no canonical OME equivalent.  Keyed by *written* channel
# position (0..C-1) to align with the Channel:0:k IDs emitted above.

.qpi_map_annotation <- function(meta, ch_meta, C) {
  if (is.null(meta)) return("")
  sl <- meta$slide %||% list()
  ii <- meta$images[[1L]]$image_info %||% list()
  cam <- ii$camera %||% list()

  kv <- list(
    description_version  = sl$description_version,
    acquisition_software = sl$acquisition_software,
    identifier           = sl$identifier,
    slide_id             = sl$slide_id,
    barcode              = sl$barcode,
    study_name           = sl$study_name,
    computer_name        = sl$computer_name,
    datetime             = sl$datetime,
    image_type           = ii$image_type,
    scan_profile_name    = ii$scan_profile_name,
    scan_mode            = ii$scan_mode,
    is_tma               = ii$is_tma,
    opal_kit_type        = ii$opal_kit_type,
    objective            = ii$objective,
    bf_lamp_type         = ii$bf_lamp_type,
    acquisition_format   = qpi_format(meta),
    camera_name          = cam$camera_name,
    camera_gain          = cam$gain,
    camera_bit_depth     = cam$bit_depth,
    channel_count        = C
  )

  # Per-channel vendor fields (no OME home).
  ch_fields <- c("is_unmixed_component", "signal_units", "objective",
                 "autofluorescence_subtracted", "responsivity",
                 "responsivity_filter_id", "responsivity_date",
                 "responsivity_filter_name", "excitation_filter_part_no",
                 "emission_filter_part_no", "bit_depth", "offset_counts",
                 "camera_orientation", "roi_x", "roi_y", "roi_width", "roi_height")
  for (k in seq_len(C)) {
    cm <- ch_meta[[k]] %||% list()
    for (f in ch_fields)
      if (!is.null(cm[[f]])) kv[[sprintf("ch%d_%s", k - 1L, f)]] <- cm[[f]]
  }

  kv <- kv[!vapply(kv, is.null, logical(1L))]
  if (length(kv) == 0L) return("")

  entries <- paste0(vapply(names(kv), function(k)
    sprintf('<M K="%s">%s</M>', .xml_escape(k),
            .xml_escape(as.character(kv[[k]]))),
    character(1L)), collapse = "")

  paste0(
    '<StructuredAnnotations>',
    '<MapAnnotation ID="Annotation:0" Namespace="qpi://vectra"><Value>',
    entries,
    '</Value></MapAnnotation></StructuredAnnotations>'
  )
}

# ---- Internal: serialise description list to XML --------------------------
#
# Output format is XML so that read_qptiff can parse <Biomarker> to recover
# channel names.  bgnorm parameters go in a <bgnorm> child element as JSON.

# Escape the three characters that are unsafe in XML text content.  Attribute
# quotes are intentionally left untouched so embedded JSON survives verbatim.
.xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# Serialise the <Biomarker>/<transform>/<bgnorm> elements shared by the
# PerkinElmerQPI blocks and the OME channel-0 metadata.  When ns_reset = TRUE
# each element declares an empty namespace (xmlns="") so it is a no-namespace
# child usable inside the OME root.
.desc_body <- function(lst, ns_reset = FALSE) {
  a         <- if (ns_reset) ' xmlns=""' else ""
  ch_name   <- lst$channel
  transform <- lst$transform
  bgnorm    <- lst$bgnorm

  bgnorm_elem <- ""
  if (!is.null(bgnorm)) {
    if (requireNamespace("jsonlite", quietly = TRUE))
      bgnorm_str <- jsonlite::toJSON(bgnorm, auto_unbox = TRUE, digits = 6L)
    else
      bgnorm_str <- paste0(
        "level=",         bgnorm$level,
        ";jsd=",          bgnorm$jsd,
        ";no_signal=",    bgnorm$no_signal,
        ";quantile_norm=",bgnorm$quantile_norm,
        if (!is.null(bgnorm$threshold)) paste0(";threshold=", bgnorm$threshold) else "",
        ";means=",  paste(round(bgnorm$parameters$means, 6L), collapse = ","),
        ";sds=",    paste(round(bgnorm$parameters$sds,   6L), collapse = ","),
        ";props=",  paste(round(bgnorm$parameters$props,  6L), collapse = ",")
      )
    bgnorm_elem <- paste0("<bgnorm", a, ">", .xml_escape(bgnorm_str), "</bgnorm>")
  }

  paste0(
    "<Biomarker", a, ">", .xml_escape(ch_name),   "</Biomarker>",
    "<transform", a, ">", .xml_escape(transform), "</transform>",
    bgnorm_elem
  )
}

.desc_to_string <- function(lst) {
  paste0("<PerkinElmerQPI>", .desc_body(lst, ns_reset = FALSE), "</PerkinElmerQPI>")
}

# ---- Internal: minimal TIFF binary writer ---------------------------------
#
# Writes a little-endian standard TIFF with one 16-bit grayscale page per
# channel.  Each page carries an ImageDescription (TIFF tag 270) populated
# with the JSON metadata built by write_qptiff().
#
# Layout per page:
#   1. Pixel data  (H x W x 2 bytes, row-major uint16 LE)
#   2. Description string (null-terminated ASCII)
#   3. IFD         (2 + 11 x 12 + 4 = 138 bytes)
#
# The file header offset points to the first page's IFD.

.write_tiff_binary <- function(path, mats, descs, H, W) {
  n_pages       <- length(mats)
  bytes_per_page <- as.integer(H) * as.integer(W) * 2L
  n_tags        <- 11L
  ifd_bytes     <- 2L + n_tags * 12L + 4L  # 138 bytes

  # Null-terminated description bytes per page (UTF-8; the OME header may hold
  # \u00b5m / \u00b5s unit symbols).
  desc_raw  <- lapply(descs, function(s) c(charToRaw(enc2utf8(s)), as.raw(0L)))
  desc_size <- vapply(desc_raw, length, integer(1L))

  # Pre-compute offsets
  pixel_off <- integer(n_pages)
  desc_off  <- integer(n_pages)
  ifd_off   <- integer(n_pages)

  off <- 8L  # header is 8 bytes
  for (p in seq_len(n_pages)) {
    pixel_off[p] <- off;  off <- off + bytes_per_page
    desc_off[p]  <- off;  off <- off + desc_size[p]
    ifd_off[p]   <- off;  off <- off + ifd_bytes
  }

  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)

  # --- TIFF header (8 bytes) ---
  writeBin(charToRaw("II"), con)                         # little-endian
  writeBin(42L,  con, size = 2L, endian = "little")      # TIFF magic
  writeBin(ifd_off[1L], con, size = 4L, endian = "little")  # offset to IFD0

  # --- Pages ---
  for (p in seq_len(n_pages)) {
    # 1. Pixel data: R matrix is column-major [H,W], TIFF wants row-major
    vals <- as.integer(as.vector(t(mats[[p]])))
    writeBin(vals, con, size = 2L, endian = "little")

    # 2. Description string (null-terminated)
    writeBin(desc_raw[[p]], con)

    # 3. IFD
    n_next <- if (p < n_pages) ifd_off[p + 1L] else 0L

    writeBin(as.integer(n_tags), con, size = 2L, endian = "little")

    .tiff_entry_long  (con, 256L, W)                     # ImageWidth
    .tiff_entry_long  (con, 257L, H)                     # ImageLength
    .tiff_entry_short (con, 258L, 16L)                   # BitsPerSample
    .tiff_entry_short (con, 259L, 1L)                    # Compression=None
    .tiff_entry_short (con, 262L, 1L)                    # PhotometricInterp=BlackIsZero
    .tiff_entry_ascii (con, 270L, desc_size[p], desc_off[p])  # ImageDescription
    .tiff_entry_long  (con, 273L, pixel_off[p])          # StripOffsets
    .tiff_entry_short (con, 277L, 1L)                    # SamplesPerPixel
    .tiff_entry_long  (con, 278L, H)                     # RowsPerStrip (all)
    .tiff_entry_long  (con, 279L, bytes_per_page)        # StripByteCounts
    .tiff_entry_short (con, 284L, 1L)                    # PlanarConfig=Chunky

    writeBin(as.integer(n_next), con, size = 4L, endian = "little")  # next IFD
  }

  invisible(NULL)
}

# IFD entry helpers (each writes exactly 12 bytes)

.tiff_entry_short <- function(con, tag, value) {
  writeBin(as.integer(tag),   con, size = 2L, endian = "little")
  writeBin(3L,                con, size = 2L, endian = "little")  # type SHORT
  writeBin(1L,                con, size = 4L, endian = "little")  # count 1
  writeBin(as.integer(value), con, size = 2L, endian = "little")
  writeBin(0L,                con, size = 2L, endian = "little")  # padding
}

.tiff_entry_long <- function(con, tag, value) {
  writeBin(as.integer(tag),   con, size = 2L, endian = "little")
  writeBin(4L,                con, size = 2L, endian = "little")  # type LONG
  writeBin(1L,                con, size = 4L, endian = "little")  # count 1
  writeBin(as.integer(value), con, size = 4L, endian = "little")
}

.tiff_entry_ascii <- function(con, tag, count, offset) {
  writeBin(as.integer(tag),    con, size = 2L, endian = "little")
  writeBin(2L,                 con, size = 2L, endian = "little")  # type ASCII
  writeBin(as.integer(count),  con, size = 4L, endian = "little")
  writeBin(as.integer(offset), con, size = 4L, endian = "little")
}
