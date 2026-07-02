## ============================================================
## QPTIFFArraySeed.R
##
## DelayedArray seed backend for Akoya QPTIFF images.
## Reads individual TIFF pages on demand, supporting:
##   - Uncompressed TIFF (compression = 1)
##   - DEFLATE / zlib compressed TIFF (compression = 8, 32946)
##   - Any other compression via tiff::readTIFF fallback
##   - Tiled and stripped TIFF layouts
##   - Standard TIFF (32-bit offsets) and BigTIFF (64-bit offsets)
##
## The seed implements the DelayedArray "seed contract":
##   dim(), dimnames(), type(), extract_array()
## ============================================================

#' @import methods
#' @importFrom DelayedArray extract_array type
NULL

# ============================================================
# QPTIFFArraySeed S4 class
# ============================================================

#' Seed class for lazy on-disk QPTIFF access
#'
#' An S4 class that satisfies the \pkg{DelayedArray} seed contract, enabling
#' lazy, block-based processing of large Akoya PhenoCycler-Fusion QPTIFF
#' images.  Each channel (TIFF page) is read from disk only when data for
#' that channel is actually requested by the \pkg{DelayedArray} framework.
#'
#' Internally the seed reads the full TIFF/BigTIFF directory chain at
#' construction time (fast: only header bytes, no pixel data) and stores
#' tile/strip layout metadata for each channel page.  Pixel data is
#' fetched per-tile (when the image is tiled) or per-strip with
#' DEFLATE/zlib decompression handled by R's built-in \code{memDecompress}.
#' Unsupported compression types fall back to \code{tiff::readTIFF}.
#'
#' Do not construct this class directly; use
#' \code{\link{read_qptiff}(path, lazy = TRUE)}.
#'
#' @slot filepath   Absolute path to the QPTIFF file.
#' @slot .dim       Integer vector \code{[H, W, C]}.
#' @slot .dimnames  List of dimnames; element 3 holds channel names.
#' @slot page_layouts List of per-channel page layout objects (one per channel)
#'   each containing the TIFF structural metadata needed to fetch that page.
#' @slot dtype      Character, \code{"integer"} (raw 16-bit values 0-65535)
#'   or \code{"double"} (normalised to [0, 1]).
#' @slot metadata   List; rich QPI metadata as returned by \code{.parse_qpi_xml}.
#' @slot level      Integer; pyramid resolution level (1 = full resolution).
#'
#' @exportClass QPTIFFArraySeed
setClass("QPTIFFArraySeed", representation(
  filepath     = "character",
  .dim         = "integer",
  .dimnames    = "list",
  page_layouts = "list",
  dtype        = "character",
  metadata     = "list",
  level        = "integer"
))

# ============================================================
# Seed contract S4 methods
# ============================================================

#' @rdname QPTIFFArraySeed-class
#' @export
setMethod("dim", "QPTIFFArraySeed", function(x) x@.dim)

#' @rdname QPTIFFArraySeed-class
#' @export
setMethod("dimnames", "QPTIFFArraySeed", function(x) {
  nms <- x@.dimnames
  if (all(vapply(nms, is.null, logical(1L)))) NULL else nms
})

#' @rdname QPTIFFArraySeed-class
#' @export
setMethod("type", "QPTIFFArraySeed", function(x) x@dtype)

#' Extract a sub-array from a QPTIFFArraySeed
#'
#' Called internally by the \pkg{DelayedArray} framework.  Reads only the
#' requested channels and the spatial region specified by \code{index} from
#' disk, keeping memory use proportional to the data actually needed.
#'
#' @param x     A \code{QPTIFFArraySeed}.
#' @param index List of length 3.  Each element is either \code{NULL} (all
#'   indices for that dimension) or a sorted integer vector of 1-based indices.
#'
#' @return An ordinary array with dimensions
#'   \code{c(length(i_rows), length(i_cols), length(i_channels))}.
#'
#' @rdname QPTIFFArraySeed-class
#' @export
setMethod("extract_array", "QPTIFFArraySeed", function(x, index) {
  h  <- x@.dim[1L]; w <- x@.dim[2L]; nc <- x@.dim[3L]

  i_rows <- index[[1L]] %||% seq_len(h)
  i_cols <- index[[2L]] %||% seq_len(w)
  i_chns <- index[[3L]] %||% seq_len(nc)

  nr <- length(i_rows); nco <- length(i_cols); nch <- length(i_chns)

  if (nr == 0L || nco == 0L || nch == 0L) {
    out <- if (x@dtype == "integer") array(NA_integer_, c(nr, nco, nch))
           else                      array(NA_real_,    c(nr, nco, nch))
    return(out)
  }

  out <- if (x@dtype == "integer") array(NA_integer_, c(nr, nco, nch))
         else                      array(NA_real_,    c(nr, nco, nch))

  for (k in seq_len(nch)) {
    layout <- x@page_layouts[[i_chns[k]]]
    pg     <- .read_qptiff_page(x@filepath, layout, i_rows, i_cols)

    if (x@dtype == "double") {
      bps     <- layout$bits_per_sample %||% 16L
      max_val <- if (bps == 8L) 255 else 65535
      pg <- pg / max_val
    }
    out[, , k] <- pg
  }

  out
})

# ============================================================
# Comprehensive IFD layout reader
# ============================================================

# TIFF tags required for reading pixel data
.PIXEL_TAGS <- c(
  "256"  = "image_width",
  "257"  = "image_height",
  "258"  = "bits_per_sample",
  "259"  = "compression",
  "273"  = "strip_offsets",
  "278"  = "rows_per_strip",
  "279"  = "strip_byte_counts",
  "322"  = "tile_width",
  "323"  = "tile_height",
  "324"  = "tile_offsets",
  "325"  = "tile_byte_counts",
  "339"  = "sample_format"
)

#' Read full page-layout metadata from every IFD in the TIFF chain
#'
#' Makes one linear pass through all TIFF/BigTIFF IFDs and collects the
#' structural tags (width, height, compression, tile/strip offsets, etc.)
#' needed to read pixel data for each page later.  No pixel data is read.
#'
#' @param path     File path.
#' @param max_pages Upper bound on number of IFDs to traverse.
#'
#' @return A list (one element per IFD) of named lists with fields:
#'   \code{image_width}, \code{image_height}, \code{bits_per_sample},
#'   \code{compression}, \code{tile_width}, \code{tile_height},
#'   \code{tile_offsets}, \code{tile_byte_counts},
#'   \code{strip_offsets}, \code{strip_byte_counts},
#'   \code{rows_per_strip}, \code{sample_format}, \code{endian}.
#' @keywords internal
.read_all_ifd_page_layouts <- function(path, max_pages = 1000L) {
  if (!file.exists(path)) return(list())
  f <- tryCatch(file(path, "rb"), error = function(e) NULL)
  if (is.null(f)) return(list())
  on.exit(close(f), add = TRUE)

  bom    <- readBin(f, "raw", n = 2L)
  endian <- if (identical(bom, as.raw(c(0x49L, 0x49L)))) "little" else "big"
  magic  <- readBin(f, "integer", n = 1L, size = 2L, signed = FALSE, endian = endian)
  bigtiff <- (magic == 43L)
  if (magic != 42L && magic != 43L) return(list())

  if (bigtiff) {
    readBin(f, "raw", n = 4L)
    ifd_offset <- .read_uint64(f, endian)
  } else {
    ifd_offset <- .ru32(f, endian)
  }

  layouts <- list()

  for (pg in seq_len(max_pages)) {
    if (is.null(ifd_offset) || ifd_offset == 0) break
    result <- tryCatch(
      .read_one_ifd_page_layout(f, ifd_offset, bigtiff, endian),
      error = function(e) NULL
    )
    if (is.null(result)) break
    layouts[[length(layouts) + 1L]] <- c(result$tags,
                                          list(endian  = endian,
                                               bigtiff = bigtiff))
    ifd_offset <- result$next_offset
  }

  layouts
}

#' Read all pixel-related tags from one IFD
#' @keywords internal
.read_one_ifd_page_layout <- function(f, offset, bigtiff, endian) {
  seek(f, offset)

  n_entries  <- if (bigtiff) .read_uint64(f, endian) else
                readBin(f, "integer", n = 1L, size = 2L, signed = FALSE, endian = endian)
  if (is.null(n_entries) || n_entries == 0L) return(NULL)

  entry_size  <- if (bigtiff) 20L else 12L
  vfield_size <- if (bigtiff)  8L else  4L
  entry_start <- seek(f)

  tag_ids  <- as.integer(names(.PIXEL_TAGS))
  tag_data <- vector("list", length(.PIXEL_TAGS))
  names(tag_data) <- names(.PIXEL_TAGS)

  for (e in seq_len(min(n_entries, 1000L))) {
    seek(f, entry_start + (e - 1L) * entry_size)
    tag <- readBin(f, "integer", n = 1L, size = 2L, signed = FALSE, endian = endian)
    if (is.null(tag) || !tag %in% tag_ids) next

    tiff_type <- readBin(f, "integer", n = 1L, size = 2L, signed = FALSE, endian = endian)
    count     <- if (bigtiff) .read_uint64(f, endian) else {
                   v <- readBin(f, "integer", n = 1L, size = 4L, signed = TRUE, endian = endian)
                   if (length(v) == 0L) next
                   if (v < 0L) as.numeric(v) + 2^32 else as.numeric(v)
                 }
    if (is.null(count)) next
    vfield <- readBin(f, "raw", n = vfield_size)

    vals <- .read_tag_values(f, tiff_type, count, vfield, bigtiff, endian)
    if (!is.null(vals))
      tag_data[[as.character(tag)]] <- vals
  }

  # Rename to friendly names
  tags <- setNames(tag_data, .PIXEL_TAGS[names(tag_data)])

  # Next IFD offset
  seek(f, entry_start + n_entries * entry_size)
  next_offset <- if (bigtiff) .read_uint64(f, endian) else .ru32(f, endian)

  list(tags = tags, next_offset = next_offset)
}

#' Read a TIFF tag's value(s) as a numeric vector
#'
#' Handles inline vs offset storage and all common TIFF types
#' (BYTE, SHORT, LONG, LONG8).  Returns a double vector for large offsets
#' to avoid signed 32-bit integer overflow.
#' @keywords internal
.read_tag_values <- function(f, tiff_type, count, vfield, bigtiff, endian) {
  vfield_size <- if (bigtiff) 8L else 4L

  type_bytes <- c(
    "1"  = 1L, "2"  = 1L, "3"  = 2L, "4"  = 4L,  "5"  = 8L,
    "6"  = 1L, "7"  = 1L, "8"  = 2L, "9"  = 4L,  "10" = 8L,
    "11" = 4L, "12" = 8L, "16" = 8L, "17" = 8L,  "18" = 8L
  )
  key <- as.character(tiff_type)
  ts  <- type_bytes[key]
  if (is.na(ts)) return(NULL)

  n_bytes <- count * ts

  if (n_bytes <= vfield_size) {
    raw_data <- vfield[seq_len(min(n_bytes, vfield_size))]
  } else {
    off <- if (bigtiff) .raw_to_uint64(vfield, endian)
           else         .raw_to_uint32(vfield, endian)
    cur <- seek(f)
    seek(f, off)
    raw_data <- readBin(f, "raw", n = n_bytes)
    seek(f, cur)
    if (length(raw_data) < n_bytes) return(NULL)
  }

  switch(key,
    "3"  = readBin(raw_data, "integer", n = count, size = 2L,
                   signed = FALSE, endian = endian),
    "4"  = {
      # uint32: read as signed, then convert negative values (> 2^31-1) to double
      vals <- readBin(raw_data, "integer", n = count, size = 4L,
                      signed = TRUE, endian = endian)
      as.numeric(vals) + ifelse(vals < 0L, 2^32, 0)
    },
    "16" = , "18" = {
      # LONG8 / IFD8 (BigTIFF): eight bytes per value
      vapply(seq_len(count), function(i) {
        .raw_to_uint64(raw_data[((i - 1L) * 8L + 1L):(i * 8L)], endian)
      }, numeric(1L))
    },
    NULL
  )
}

# ============================================================
# Page reading dispatcher
# ============================================================

#' Read one TIFF page as a numeric matrix
#'
#' Routes to the appropriate reader based on the page layout: tiled DEFLATE,
#' tiled uncompressed, stripped DEFLATE, stripped uncompressed, or a fallback
#' through \code{tiff::readTIFF}.
#'
#' @param path   File path.
#' @param layout Page layout list from \code{.read_all_ifd_page_layouts}.
#' @param rows   Integer vector of 1-based row indices to return, or \code{NULL}
#'   for all rows.
#' @param cols   Integer vector of 1-based column indices to return, or
#'   \code{NULL} for all columns.
#'
#' @return Integer matrix \code{[length(rows), length(cols)]}.
#' @keywords internal
.read_qptiff_page <- function(path, layout, rows = NULL, cols = NULL) {
  h   <- as.integer(layout$image_height %||% 0L)
  w   <- as.integer(layout$image_width  %||% 0L)
  if (h == 0L || w == 0L)
    stop("Cannot determine page dimensions from TIFF layout.")

  if (is.null(rows)) rows <- seq_len(h)
  if (is.null(cols)) cols <- seq_len(w)

  bps         <- as.integer(layout$bits_per_sample %||% 16L)
  compression <- as.integer(layout$compression %||% 1L)

  # Dispatch by layout type; tiled/stripped readers handle all compressions
  # (native DEFLATE; everything else via per-tile minimal-TIFF decode)
  if (!is.null(layout$tile_offsets))
    return(.read_tiled_page(path, layout, h, w, bps, compression, rows, cols))
  if (!is.null(layout$strip_offsets))
    return(.read_stripped_page(path, layout, h, w, bps, compression, rows, cols))

  # Fallback: no tile/strip layout available - use tiff package
  .read_page_fallback(path, layout, rows, cols)
}

# ============================================================
# Tiled page reader
# ============================================================

#' Read a tiled TIFF page (handles uncompressed + DEFLATE)
#' @keywords internal
.read_tiled_page <- function(path, layout, h, w, bps, compression, rows, cols) {
  tile_w    <- as.integer(layout$tile_width)
  tile_h    <- as.integer(layout$tile_height)
  tile_offs <- as.numeric(layout$tile_offsets)
  tile_bc   <- as.numeric(layout$tile_byte_counts)
  endian    <- layout$endian %||% "little"

  n_tcols   <- ceiling(w / tile_w)
  n_trows   <- ceiling(h / tile_h)

  out <- matrix(0L, nrow = length(rows), ncol = length(cols))

  f <- tryCatch(file(path, "rb"), error = function(e) NULL)
  if (is.null(f)) stop("Cannot open file: ", path)
  on.exit(close(f), add = TRUE)

  bpp <- bps %/% 8L       # bytes per pixel
  n_tile_bytes <- tile_h * tile_w * bpp   # expected decompressed bytes per tile

  for (tr in seq_len(n_trows)) {
    # Row range covered by this tile row (1-based)
    tile_r0 <- (tr - 1L) * tile_h + 1L
    tile_r1 <- min(tr * tile_h, h)

    rows_in <- rows[rows >= tile_r0 & rows <= tile_r1]
    if (length(rows_in) == 0L) next
    out_rows <- match(rows_in, rows)
    local_rows <- rows_in - tile_r0 + 1L

    for (tc in seq_len(n_tcols)) {
      tile_c0 <- (tc - 1L) * tile_w + 1L
      tile_c1 <- min(tc * tile_w, w)

      cols_in <- cols[cols >= tile_c0 & cols <= tile_c1]
      if (length(cols_in) == 0L) next
      out_cols <- match(cols_in, cols)
      local_cols <- cols_in - tile_c0 + 1L

      tidx <- (tr - 1L) * n_tcols + tc
      if (tidx > length(tile_offs)) next

      off    <- tile_offs[tidx]
      nbytes <- as.integer(tile_bc[tidx])

      seek(f, off)
      compressed <- readBin(f, "raw", n = nbytes)

      raw_pixels <- .decompress_tiff_data(compressed, compression, n_tile_bytes)
      if (is.null(raw_pixels)) {
        tile_mat <- .decode_tile_via_tiff(compressed, tile_w, tile_h, bps,
                                          compression, endian)
        if (is.null(tile_mat)) next
      } else {
        tile_mat <- .raw_to_int_matrix(raw_pixels, tile_h, tile_w, bps, endian)
      }

      out[out_rows, out_cols] <- tile_mat[local_rows, local_cols, drop = FALSE]
    }
  }

  out
}

# ============================================================
# Stripped page reader
# ============================================================

#' Read a stripped TIFF page (handles uncompressed + DEFLATE)
#' @keywords internal
.read_stripped_page <- function(path, layout, h, w, bps, compression, rows, cols) {
  strip_offs  <- as.numeric(layout$strip_offsets)
  strip_bc    <- as.numeric(layout$strip_byte_counts)
  rps         <- as.integer(layout$rows_per_strip %||% h)
  endian      <- layout$endian %||% "little"
  n_strips    <- length(strip_offs)
  bpp         <- bps %/% 8L

  out <- matrix(0L, nrow = length(rows), ncol = length(cols))

  f <- tryCatch(file(path, "rb"), error = function(e) NULL)
  if (is.null(f)) stop("Cannot open file: ", path)
  on.exit(close(f), add = TRUE)

  for (s in seq_len(n_strips)) {
    strip_r0 <- (s - 1L) * rps + 1L
    strip_r1 <- min(s * rps, h)

    rows_in <- rows[rows >= strip_r0 & rows <= strip_r1]
    if (length(rows_in) == 0L) next

    n_strip_rows  <- strip_r1 - strip_r0 + 1L
    n_strip_bytes <- n_strip_rows * w * bpp

    seek(f, strip_offs[s])
    compressed <- readBin(f, "raw", n = as.integer(strip_bc[s]))

    raw_pixels <- .decompress_tiff_data(compressed, compression, n_strip_bytes)
    if (is.null(raw_pixels)) {
      strip_mat <- .decode_tile_via_tiff(compressed, w, n_strip_rows, bps,
                                         compression, endian)
      if (is.null(strip_mat)) next
    } else {
      strip_mat <- .raw_to_int_matrix(raw_pixels, n_strip_rows, w, bps, endian)
    }

    out_rows   <- match(rows_in, rows)
    local_rows <- rows_in - strip_r0 + 1L
    out[out_rows, ] <- strip_mat[local_rows, cols, drop = FALSE]
  }

  out
}

# ============================================================
# Fallback reader
# ============================================================

#' Fallback page reader via tiff::readTIFF
#'
#' Used when compression is not natively supported (LZW, JPEG, etc.).
#' Loads ALL pages into memory then extracts the requested page; this
#' is memory-intensive for large multi-channel files.
#'
#' @keywords internal
.read_page_fallback <- function(path, layout, rows, cols) {
  pg_idx <- layout$page_idx %||% 1L

  if (requireNamespace("ijtiff", quietly = TRUE)) {
    result <- tryCatch({
      img <- ijtiff::read_tif(path, frames = pg_idx, msg = FALSE)
      # ijtiff returns [H, W, C, F]; take first channel and frame
      m <- img[, , 1L, 1L, drop = TRUE]
      m[rows, cols, drop = FALSE]
    }, error = function(e) NULL)
    if (!is.null(result)) return(result)
  }

  warning("Loading all QPTIFF pages to read page ", pg_idx,
          " (unsupported compression - install 'ijtiff' for efficient single-page reads).")
  all_pages <- tiff::readTIFF(path, all = TRUE, as.is = TRUE)
  if (pg_idx > length(all_pages))
    stop("Page index ", pg_idx, " exceeds page count ", length(all_pages))
  pg <- all_pages[[pg_idx]]
  if (length(dim(pg)) == 3L) pg <- pg[, , 1L]
  pg[rows, cols, drop = FALSE]
}

# ============================================================
# Decompression and byte conversion helpers
# ============================================================

#' Decompress one TIFF tile or strip
#'
#' Handles uncompressed (compression = 1) and DEFLATE/zlib (8, 32946).
#' For DEFLATE, uses \code{memDecompress(type="gzip")} which correctly
#' handles zlib-wrapped deflate streams (0x78 0x9C / 0x78 0xDA headers)
#' produced by the Akoya PhenoCycler-Fusion software.
#'
#' @param raw_bytes    Raw vector of compressed bytes from file.
#' @param compression  TIFF compression code.
#' @param n_bytes_out  Expected number of bytes after decompression.
#'
#' @return Raw vector of decompressed bytes.
#' @keywords internal
.decompress_tiff_data <- function(raw_bytes, compression, n_bytes_out) {
  switch(as.character(compression),
    "1"     = raw_bytes,
    "8"     = ,  # Deflate (zlib)
    "32946" = {  # Adobe Deflate (also zlib)
      out <- tryCatch(memDecompress(raw_bytes, type = "gzip"),
                      error = function(e) stop("DEFLATE decompression failed: ",
                                               conditionMessage(e)))
      if (length(out) < n_bytes_out)
        stop("Decompressed ", length(out), " bytes; expected ", n_bytes_out)
      out
    },
    NULL  # Return NULL for unsupported compression; caller uses .decode_tile_via_tiff
  )
}

#' Decode a single compressed tile or strip via an in-memory TIFF wrapper
#'
#' Wraps the raw compressed bytes in a minimal single-strip TIFF structure and
#' decodes via \code{tiff::readTIFF()}.  This delegates all decompression
#' (LZW, JPEG, PackBits, etc.) to libtiff without loading the entire source file.
#'
#' @param raw_bytes   Raw vector of compressed tile/strip bytes.
#' @param w,h         Tile/strip width and height in pixels.
#' @param bps         Bits per sample (8, 16, or 32).
#' @param compression TIFF compression code.
#' @param endian      Byte order for the in-memory TIFF header.
#'
#' @return Integer matrix \code{[h, w]}, or \code{NULL} on error.
#' @keywords internal
.decode_tile_via_tiff <- function(raw_bytes, w, h, bps, compression,
                                  endian = "little") {
  n_tags   <- 8L
  data_off <- 8L + 2L + n_tags * 12L + 4L  # TIFF header + IFD = 110 bytes

  to2 <- function(x) writeBin(as.integer(x), raw(), size = 2L, endian = endian)
  to4 <- function(x) writeBin(as.integer(x), raw(), size = 4L, endian = endian)
  ifd_entry <- function(tag, type, cnt, val)
    c(to2(tag), to2(type), to4(cnt), to4(val))

  hdr <- c(as.raw(c(0x49, 0x49, 0x2A, 0x00)), to4(8L))  # "II" + 42 + IFD@8

  # IFD tags must be in ascending numeric order
  ifd <- c(
    to2(n_tags),
    ifd_entry(256L, 4L, 1L, w),                    # ImageWidth  (LONG)
    ifd_entry(257L, 4L, 1L, h),                    # ImageLength (LONG)
    ifd_entry(258L, 3L, 1L, bps),                  # BitsPerSample (SHORT)
    ifd_entry(259L, 3L, 1L, compression),          # Compression  (SHORT)
    ifd_entry(262L, 3L, 1L, 1L),                   # PhotometricInterpretation (SHORT)
    ifd_entry(273L, 4L, 1L, data_off),             # StripOffsets (LONG)
    ifd_entry(278L, 4L, 1L, h),                    # RowsPerStrip (LONG)
    ifd_entry(279L, 4L, 1L, length(raw_bytes)),    # StripByteCounts (LONG)
    to4(0L)                                         # next IFD = 0
  )

  tiff_data <- c(hdr, ifd, raw_bytes)

  mat <- tryCatch(
    suppressWarnings(tiff::readTIFF(tiff_data, all = FALSE, as.is = TRUE)),
    error = function(e) NULL
  )
  if (is.null(mat)) return(NULL)
  if (length(dim(mat)) == 3L) mat <- mat[, , 1L]           # drop colour dim
  # tiff may normalise 16-bit to [0,1] for some compressors - restore raw range
  if (is.double(mat) && bps == 16L && max(mat, na.rm = TRUE) <= 1.0)
    mat <- round(mat * 65535)
  if (!is.matrix(mat)) return(NULL)
  mat
}

#' Convert raw pixel bytes to an integer matrix
#'
#' @param raw_bytes Raw vector (decompressed pixel data, row-major).
#' @param nr        Number of rows.
#' @param nc        Number of columns.
#' @param bps       Bits per sample (8 or 16).
#' @param endian    Byte order ("little" or "big").
#'
#' @return Integer matrix \code{[nr, nc]}.
#' @keywords internal
.raw_to_int_matrix <- function(raw_bytes, nr, nc, bps, endian) {
  n <- nr * nc
  if (bps == 16L) {
    vals <- readBin(raw_bytes[seq_len(n * 2L)], "integer", n = n,
                    size = 2L, signed = FALSE, endian = endian)
  } else if (bps == 8L) {
    vals <- as.integer(raw_bytes[seq_len(n)])
  } else if (bps == 32L) {
    vals <- readBin(raw_bytes[seq_len(n * 4L)], "integer", n = n,
                    size = 4L, signed = TRUE, endian = endian)
    vals <- as.numeric(vals) + ifelse(vals < 0L, 2^32, 0)
  } else {
    stop("Unsupported bits-per-sample: ", bps)
  }
  matrix(vals, nrow = nr, ncol = nc, byrow = TRUE)
}
