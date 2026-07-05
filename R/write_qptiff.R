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
  descs[1L] <- .build_ome_description(chs, H, W, C, desc_lists[[1L]])

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

.build_ome_description <- function(chs, H, W, C, ch0_desc_list) {
  ns <- "http://www.openmicroscopy.org/Schemas/OME/2016-06"

  channels <- paste0(vapply(seq_len(C), function(k)
    sprintf('<Channel ID="Channel:0:%d" Name="%s" SamplesPerPixel="1"/>',
            k - 1L, .xml_escape(chs[k])),
    character(1L)), collapse = "")

  tiffdata <- paste0(vapply(seq_len(C), function(k)
    sprintf('<TiffData FirstC="%d" FirstZ="0" FirstT="0" IFD="%d" PlaneCount="1"/>',
            k - 1L, k - 1L),
    character(1L)), collapse = "")

  paste0(
    '<?xml version="1.0" encoding="UTF-8"?>',
    '<OME xmlns="', ns, '" ',
    'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" ',
    'xsi:schemaLocation="', ns, ' ', ns, '/ome.xsd">',
    '<Image ID="Image:0" Name="bgnormR">',
    '<Pixels ID="Pixels:0" DimensionOrder="XYCZT" Type="uint16" ',
    'SizeX="', W, '" SizeY="', H, '" SizeC="', C, '" SizeZ="1" SizeT="1" ',
    'Interleaved="false" BigEndian="false" SignificantBits="16">',
    channels, tiffdata,
    '</Pixels></Image>',
    .desc_body(ch0_desc_list, ns_reset = TRUE),
    '</OME>'
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

  # Null-terminated description bytes per page
  desc_raw  <- lapply(descs, function(s) c(charToRaw(s), as.raw(0L)))
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
