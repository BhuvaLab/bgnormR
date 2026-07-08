## ============================================================
## read_qptiff.R
##
## QPTIFF reader ported from Python bioio-tifffile
## (rtubelleza/bioio-tifffile, feature/read-qptiffs-rich-reduce-ome).
##
## Handles three live PerkinElmer / Akoya format variants:
##   brightfield     - H&E / IHC; RGB pages
##   polaris_scanband - Vectra / Polaris / OPAL; channel metadata in
##                     <ScanBands-i> XML elements
##   fusion_paged    - Akoya PhenoCycler-Fusion; per-page XML with JSON
##                     <ScanProfile>
##
## No Java dependency; uses tiff + xml2 (both in Suggests).
## ============================================================

# ----------- Format constants (internal) ----------------------------------

.QPTIFF_BF  <- "brightfield"
.QPTIFF_PSB <- "polaris_scanband"
.QPTIFF_FP  <- "fusion_paged"
.QPTIFF_UNK <- "unknown"

.BIOMARKER_TAGS   <- c("Biomarker", "BioMarker", "StainName", "Marker", "Name")
.FLUOROPHORE_TAGS <- c("Fluorophore", "Fluor")

# TIFF tag numbers
.TAG_IMAGE_DESC    <- 270L
.TAG_SAMPLES_PIX   <- 277L

# ============================================================
# Public API
# ============================================================

#' Read a multiplex image (QPTIFF, OME-TIFF, or OME-Zarr)
#'
#' Reads a multiplex image into R without Java, auto-detecting the container:
#' \itemize{
#'   \item \strong{QPTIFF} from Akoya PhenoCycler-Fusion (formerly CODEX),
#'     Cell DIVE, or Vectra / Polaris scanners - brightfield RGB, Polaris
#'     ScanBand XML, and Fusion per-page JSON+XML variants.
#'   \item \strong{OME-TIFF} (as written by tifffile, Bio-Formats, or QuPath) -
#'     channel names are read from the \code{<Channel Name="...">} attributes of
#'     the OME-XML.
#'   \item \strong{OME-Zarr} (OME-NGFF) stores - requires the \pkg{Rarr}
#'     package; channel names / colours are read from \code{omero.channels}.
#' }
#'
#' Channel names and rich per-channel metadata (fluorophore, exposure time,
#' wavelengths, filters, colour) are parsed into an OME-organised
#' \code{\link{QPTIFFMetadata}} object, available via \code{\link{metadata}}.
#'
#' The native, Java-free TIFF/QPTIFF reader and writer implemented here was
#' translated from the
#' \href{https://github.com/rtubelleza/bioio-tifffile/tree/feature/read-qptiffs-rich}{bioio-tifffile fork}
#' by Rafael Tubelleza.
#'
#' @param path       Character, path to the image.  A QPTIFF / OME-TIFF file, or
#'   an OME-Zarr store (a directory, or a path ending in \code{.zarr}).
#' @param channels   Character vector of channel names to load, an integer
#'   vector of 1-based channel indices, or \code{NULL} (default) to load all
#'   channels.
#' @param level      Integer, pyramid resolution level.  \code{1} = full
#'   resolution (default), \code{2} = half resolution, etc.  OME-TIFF
#'   sub-resolution pyramids (SubIFDs) are not supported; use \code{level = 1}.
#' @param as_integer Logical; return raw 16-bit integers (0-65535) rather than
#'   normalised \code{[0, 1]} doubles?  Default \code{TRUE}.
#' @param lazy       Logical; if \code{TRUE} return a
#'   \code{\link[DelayedArray]{DelayedArray}} backed by a
#'   \code{\linkS4class{QPTIFFArraySeed}} that reads pages from disk on demand.
#'   If \code{FALSE} (default) load all requested channels into memory and
#'   return a \code{\link{QPTIFFImage}}.
#'
#' @return
#' \describe{
#'   \item{\code{lazy = FALSE}}{A \code{QPTIFFImage} - a 3-D numeric array
#'     \code{[height, width, channels]} with class \code{c("QPTIFFImage",
#'     "array")}.  Channel names are stored in \code{dimnames(img)[[3]]}.
#'     Standard array subscripting works: \code{img[, , "DAPI"]} extracts a
#'     single channel as a matrix; \code{img[1:512, 1:512, ]} crops spatially.
#'     Rich metadata is in \code{attr(img, "metadata")}.}
#'   \item{\code{lazy = TRUE}}{A \code{DelayedArray} with
#'     \code{dim = c(H, W, C)}.  Individual channel pages are read from disk
#'     only when accessed.  The seed is a \code{QPTIFFArraySeed} accessible via
#'     \code{DelayedArray::seed(arr)}.}
#' }
#'
#' @importFrom utils head
#' @importFrom DelayedArray DelayedArray
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' dim(img)                      # c(H, W, n_channels)
#' names(img)                    # channel names
#' cd20 <- img[, , "CD20"]      # extract one channel as a 2-D matrix
#'
#' # Load specific channels only
#' img2 <- read_qptiff(path, channels = c("CD20", "PanCK"))
#' names(img2)
#'
#' \donttest{
#' # Lazy / out-of-core load (backed by DelayedArray)
#' arr <- read_qptiff(path, lazy = TRUE)
#' }
read_qptiff <- function(path, channels = NULL, level = 1L,
                         as_integer = TRUE, lazy = FALSE) {
  if (!is.character(path) || length(path) != 1L)
    stop("'path' must be a single character string.")
  if (!file.exists(path))
    stop("File not found: ", path)

  level <- as.integer(level)
  stopifnot(level >= 1L)
  stopifnot(is.logical(lazy), length(lazy) == 1L)

  # --- 0. Dispatch: OME-Zarr (NGFF) store --------------------------------
  # An OME-Zarr store is a directory (or a .zarr path); route it to the
  # zarr reader before any TIFF-specific handling.
  if (.is_ome_zarr_store(path))
    return(.read_ome_zarr(path, channels = channels, level = level,
                          as_integer = as_integer, lazy = lazy))

  if (!requireNamespace("tiff", quietly = TRUE))
    stop("Package 'tiff' is required.  Install with: install.packages('tiff')")
  if (!requireNamespace("xml2", quietly = TRUE))
    stop("Package 'xml2' is required.  Install with: install.packages('xml2')")

  # --- 1. Read all IFD ImageDescription XMLs + SamplesPerPixel ----------
  message("Reading TIFF directory structure ...")
  ifd_info <- .read_all_ifd_info(path)

  if (length(ifd_info$descriptions) == 0L)
    stop("No TIFF pages found in: ", path)

  # --- 1b. Dispatch: standard OME-TIFF (OME-XML root) --------------------
  # QuPath-style OME-TIFFs carry an <OME> root in tag 270 with channel names
  # as <Channel Name="..."> attributes; the QPI parser cannot read those, so
  # route to the dedicated OME reader.
  if (.is_ome_xml(ifd_info$descriptions[[1L]]))
    return(.read_ome_tiff(path, ifd_info, channels = channels, level = level,
                          as_integer = as_integer, lazy = lazy))

  # --- 2. Determine format, channel count, pyramid levels ---------------
  layout <- .detect_layout(ifd_info$descriptions, ifd_info$samples_per_pixel)
  n_ch    <- layout$n_channels
  n_lev   <- layout$n_levels
  is_rgb  <- layout$is_rgb

  if (n_ch == 0L)
    stop("Could not determine channel count from QPTIFF metadata.")

  # --- 3. Parse rich metadata from first-page XML + per-page XMLs -------
  first_xml <- ifd_info$descriptions[[1L]]
  pp_xmls   <- ifd_info$descriptions[seq_len(n_ch)]

  parsed <- .parse_qpi_xml(
    xml_str       = first_xml,
    per_page_xmls = pp_xmls,
    n_channels    = n_ch,
    is_rgb        = is_rgb
  )

  # Promote format to "brightfield" for is_rgb files whose XML didn't parse
  fmt <- parsed$format
  if (is_rgb && fmt == .QPTIFF_UNK) fmt <- .QPTIFF_BF

  # --- 3b. Separate protein name from fluorescence dye suffix ---------------
  # Compound biomarker names like "CD3e-AF647" are split into a protein name
  # ("CD3e") and a dye name ("AF647").  ch$name is updated to the protein
  # name; the extracted dye is stored in ch$dye_from_name.  Brightfield
  # channels (R/G/B) are left unchanged.
  channels_meta <- parsed$channels
  if (!is_rgb) {
    channels_meta <- lapply(channels_meta, function(ch) {
      split            <- .split_protein_dye(ch$name, ch$fluorophore)
      ch$name          <- split$protein
      ch$dye_from_name <- split$dye
      ch
    })
  }

  all_channel_names <- vapply(channels_meta,
                              function(ch) ch$name %||% NA_character_,
                              character(1L))

  # Channel-name fallbacks.  Too few parsed entries -> synthesize the whole set;
  # otherwise keep the parsed names and fill only the channels that lacked one.
  if (length(all_channel_names) < n_ch) {
    all_channel_names <- if (is_rgb && n_ch == 3L) c("Red", "Green", "Blue")
      else paste0(if (is_rgb) "Sample_" else "Channel_", seq_len(n_ch) - 1L)
  } else if (anyNA(all_channel_names)) {
    miss <- which(is.na(all_channel_names))
    all_channel_names[miss] <- paste0(if (is_rgb) "Sample_" else "Channel_",
                                       miss - 1L)
  }

  # Normalise the channel-metadata list to n_ch entries, aligned to the names.
  channels_meta <- lapply(seq_len(n_ch), function(i) {
    ch <- if (i <= length(channels_meta)) channels_meta[[i]] else list()
    ch$index          <- ch$index %||% (i - 1L)
    ch$name           <- all_channel_names[i]
    ch$is_brightfield <- ch$is_brightfield %||% is_rgb
    ch
  })

  # --- 4. Select channels -----------------------------------------------
  sel      <- .select_channels(channels, all_channel_names)
  ch_idx   <- sel$ch_idx
  ch_names <- sel$ch_names

  # --- 4b. Assemble hierarchical metadata for the *loaded* channels ------
  base_px <- parsed$image_info$scan_resolution$base_pixel_size_um
  meta <- .single_scene_metadata(
    slide      = parsed$slide,
    image_info = parsed$image_info,
    channels   = channels_meta[ch_idx],
    scales     = .build_scales(n_lev, base_px),
    format     = fmt,
    raw_xml    = parsed$raw_xml
  )

  # --- 5. Read all IFD page layouts (fast: no pixel data) ---------------
  message("Reading IFD page layouts ...")
  all_layouts_raw <- .read_all_ifd_page_layouts(path)

  # --- 6. Compute channel -> absolute page index mapping -----------------
  if (is_rgb) {
    # Brightfield: all channels are samples within one page per pyramid level.
    # page_indices[k] = level for all k; the eager/lazy path splits samples.
    page_indices <- rep(level, length(ch_idx))
  } else {
    page_indices <- tryCatch(
      .compute_data_page_indices(all_layouts_raw, n_ch, level, ch_idx),
      error = function(e) {
        # Fallback: arithmetic indexing (correct when no thumbnail pages)
        (level - 1L) * n_ch + ch_idx
      }
    )
  }

  # --- 6.5  Validate that the requested level's data is in the file --------
  if (!is_rgb && length(page_indices) > 0L) {
    file_sz <- file.size(path)
    .fdoff <- function(l) {
      if (!is.null(l$tile_offsets)  && length(l$tile_offsets)  > 0L)
        return(as.numeric(l$tile_offsets[[1L]]))
      if (!is.null(l$strip_offsets) && length(l$strip_offsets) > 0L)
        return(as.numeric(l$strip_offsets[[1L]]))
      NA_real_
    }
    first_off <- .fdoff(all_layouts_raw[[page_indices[1L]]])
    if (!is.na(first_off) && first_off >= file_sz) {
      avail <- Filter(function(lv) {
        pg <- tryCatch(
          .compute_data_page_indices(all_layouts_raw, n_ch, lv, 1L),
          error = function(e) integer(0L)
        )
        if (length(pg) == 0L) return(FALSE)
        off <- .fdoff(all_layouts_raw[[pg[1L]]])
        !is.na(off) && off < file_sz
      }, seq_len(n_lev))
      stop(
        "Level ", level, " pixel data is not accessible in this file ",
        "(file may be truncated).\n",
        "  File size:       ", file_sz, " bytes (", round(file_sz / 1e9, 2), " GB)\n",
        "  Level ", level, " data offset: ", round(first_off / 1e6, 1), " MB\n",
        if (length(avail) > 0L)
          paste0("  Levels with accessible pixel data: ",
                 paste(avail, collapse = ", "), "\n",
                 "  Hint: try level = ", max(avail))
        else
          "  No levels with accessible pixel data found in this file."
      )
    }
  }

  # --- 7. Branch: lazy (DelayedArray) vs eager (QPTIFFImage) -----------
  if (lazy) {
    .read_qptiff_lazy(path, ch_idx, ch_names, page_indices,
                      all_layouts_raw, level, meta, as_integer)
  } else {
    .read_qptiff_eager(path, ch_idx, ch_names, page_indices,
                       all_layouts_raw, meta, as_integer, is_rgb)
  }
}

# ---- private: eager path -----------------------------------------------

.read_qptiff_eager <- function(path, ch_idx, ch_names, page_indices,
                                all_layouts, meta, as_integer,
                                is_rgb = FALSE) {
  message("Loading ", length(ch_idx), " channel(s) ...")

  if (is_rgb) {
    # Brightfield: all channels (R/G/B) are interleaved samples in one page.
    rgb_page_idx <- page_indices[1L]

    # Guard: JPEG / LZW compressed RGB pages cannot be decoded by our reader
    # and crash tiff::readTIFF on some platforms (Apple Silicon bus error).
    rgb_layout <- all_layouts[[rgb_page_idx]]
    rgb_comp   <- as.integer(rgb_layout$compression %||% 1L)
    if (!rgb_comp %in% c(1L, 8L, 32946L)) {
      rgb_h <- as.integer(rgb_layout$image_height %||% 0L)
      rgb_w <- as.integer(rgb_layout$image_width  %||% 0L)
      stop("JPEG/LZW compressed RGB TIFF (compression=", rgb_comp,
           ") cannot be decoded in eager mode.\n",
           "  Image dimensions: ", rgb_h, " x ", rgb_w,
           " x ", length(ch_idx), " channels.\n",
           "  Use lazy = TRUE for metadata access without pixel loading.")
    }

    all_pgs <- tryCatch(
      tiff::readTIFF(path, all = TRUE, convert = FALSE,
                     as.is = as_integer, native = FALSE),
      error = function(e) stop("tiff::readTIFF failed: ", conditionMessage(e))
    )
    if (!is.list(all_pgs)) all_pgs <- list(all_pgs)
    if (rgb_page_idx > length(all_pgs))
      stop("Level page (", rgb_page_idx, ") not found in file.")
    pg <- all_pgs[[rgb_page_idx]]
    if (length(dim(pg)) == 3L) {
      h <- dim(pg)[1L]; w <- dim(pg)[2L]
      arr <- array(NA_real_, dim = c(h, w, length(ch_idx)))
      for (k in seq_along(ch_idx)) arr[, , k] <- pg[, , ch_idx[k]]
    } else {
      h <- nrow(pg); w <- ncol(pg)
      arr <- array(NA_real_, dim = c(h, w, length(ch_idx)))
      for (k in seq_along(ch_idx)) arr[, , k] <- pg
    }
    dimnames(arr) <- list(NULL, NULL, ch_names)
    return(.new_QPTIFFImage(arr, meta))
  }

  # Fluorescence: each channel is a separate grayscale TIFF page.
  # Load pages individually to avoid reading the entire file into memory.
  fl  <- all_layouts[[page_indices[1L]]]
  h   <- as.integer(fl$image_height %||% 0L)
  w   <- as.integer(fl$image_width  %||% 0L)
  if (h == 0L || w == 0L) {
    pg0 <- tryCatch(tiff::readTIFF(path, all = FALSE, as.is = TRUE),
                    error = function(e) NULL)
    if (!is.null(pg0)) { h <- nrow(pg0); w <- ncol(pg0) }
  }

  bps     <- as.integer(fl$bits_per_sample %||% 16L)
  max_val <- if (bps == 8L) 255.0 else 65535.0

  arr <- array(NA_real_, dim = c(h, w, length(ch_idx)))
  for (k in seq_along(ch_idx)) {
    layout <- all_layouts[[page_indices[k]]]
    pg <- tryCatch(
      .read_qptiff_page(path, layout),
      error = function(e) {
        warning("Page ", page_indices[k], ": ", conditionMessage(e),
                "\nReturning NA for channel '", ch_names[k], "'.")
        matrix(NA_real_, nrow = h, ncol = w)
      }
    )
    if (!as_integer) pg <- pg / max_val
    arr[, , k] <- pg
  }
  dimnames(arr) <- list(NULL, NULL, ch_names)
  .new_QPTIFFImage(arr, meta)
}

# ---- private: lazy path ------------------------------------------------

.read_qptiff_lazy <- function(path, ch_idx, ch_names, page_indices,
                               all_layouts, level, meta, as_integer) {
  if (length(all_layouts) > 0L && max(page_indices) > length(all_layouts))
    stop("Requested level ", level, " exceeds available pages (",
         length(all_layouts), ").")

  selected_layouts <- lapply(seq_along(ch_idx), function(k) {
    l <- all_layouts[[page_indices[k]]]
    l$page_idx <- page_indices[k]
    l
  })

  fl <- selected_layouts[[1L]]
  h  <- as.integer(fl$image_height %||% 0L)
  w  <- as.integer(fl$image_width  %||% 0L)

  if (h == 0L || w == 0L) {
    message("Reading first page to determine dimensions ...")
    pg <- tiff::readTIFF(path, all = FALSE, as.is = TRUE)
    h  <- nrow(pg); w <- ncol(pg)
  }

  seed <- new("QPTIFFArraySeed",
    filepath     = normalizePath(path),
    .dim         = c(h, w, length(ch_idx)),
    .dimnames    = list(NULL, NULL, ch_names),
    page_layouts = selected_layouts,
    dtype        = if (as_integer) "integer" else "double",
    metadata     = meta,
    level        = level
  )

  da <- DelayedArray::DelayedArray(seed)
  .new_QPTIFFImage_lazy(da, metadata = meta)
}

# ============================================================
# QPTIFFImage S3 class  (3-D array [H x W x C] with named channels)
# ============================================================

#' QPTIFFImage: an in-memory or on-disk multi-channel image
#'
#' An S3 class representing a \code{[H x W x C]} multiplex image, either
#' eagerly loaded as a plain 3-D array or lazily backed by a
#' \code{DelayedArray} (see \code{\link{QPTIFFArraySeed-class}}). Created by
#' \code{\link{read_qptiff}} or \code{\link{as.QPTIFFImage}}, and consumed by
#' \code{\link{bgnorm_pixels}}, \code{\link{plot_qptiff}}, and
#' \code{\link{write_qptiff}}.
#'
#' @return Not applicable; documents the \code{QPTIFFImage} class structure.
#' @name QPTIFFImage
NULL

#' Eager constructor for QPTIFFImage (in-memory array backing)
#'
#' @param arr            Numeric array of dimension \code{c(H, W, C)}.
#' @param metadata       List of slide / image / channel metadata.
#' @param bgnorm_results Named list of \code{BgnormResult} objects (one per
#'   channel), or \code{NULL}.
#' @keywords internal
.new_QPTIFFImage <- function(arr, metadata = list(), bgnorm_results = NULL) {
  structure(arr, class = c("QPTIFFImage", "array"),
            metadata = metadata, bgnorm_results = bgnorm_results)
}

#' Lazy constructor for QPTIFFImage (DelayedArray backing)
#'
#' @param da             A \code{DelayedArray} from a \code{QPTIFFArraySeed}.
#' @param metadata       List of slide / image / channel metadata.
#' @param bgnorm_results Named list of \code{BgnormResult} objects, or \code{NULL}.
#' @keywords internal
.new_QPTIFFImage_lazy <- function(da, metadata = list(), bgnorm_results = NULL) {
  structure(list(.da = da), class = "QPTIFFImage",
            metadata = metadata, bgnorm_results = bgnorm_results)
}

#' True when a QPTIFFImage is backed by a DelayedArray (lazy / on-disk)
#' @keywords internal
.is_lazy_qptiff <- function(x) is.list(x) && !is.null(x$.da)


#' @export
print.QPTIFFImage <- function(x, ...) {
  d   <- dim(x)
  chs <- dimnames(x)[[3L]]
  m   <- attr(x, "metadata")
  br  <- attr(x, "bgnorm_results")
  cat("QPTIFFImage\n")
  if (.is_lazy_qptiff(x)) cat("  Mode      : lazy (on-disk)\n")
  cat("  Dimensions:", d[1L], "x", d[2L], "(H x W)\n")
  cat("  Channels  :", d[3L], "\n")
  if (!is.null(chs))
    cat("  Names     :", paste(head(chs, 6L), collapse = ", "),
        if (length(chs) > 6L) "..." else "", "\n")
  if (!is.null(m)) {
    cat("  Format    :", qpi_format(m) %||% "unknown", "\n")
    px <- qpi_pixel_size_um(m)
    if (!is.null(px))
      cat("  Pixel size:", px, "um\n")
    cat("  Levels    :", qpi_n_levels(m), "\n")
  }
  if (!is.null(br))
    cat("  bgnorm    : yes (", length(br), " channel(s))\n", sep = "")
  invisible(x)
}

#' Per-channel bgnorm results stored in a QPTIFFImage
#'
#' Returns the named list of \code{\link{BgnormResult}} objects attached to a
#' \code{\link{QPTIFFImage}} by \code{\link{bgnorm_pixels}}, or \code{NULL} if
#' the image has not been background-normalised.
#'
#' @param x A \code{QPTIFFImage}.
#' @param ... Unused.
#' @return A named list of \code{BgnormResult} objects (one per processed
#'   channel), or \code{NULL}.
#' @seealso \code{\link{bgnorm_pixels}}
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' img <- read_qptiff(path)
#' adj <- bgnorm_pixels(img, markers = names(img)[1])
#' bgnorm_results(adj)
bgnorm_results <- function(x, ...) UseMethod("bgnorm_results")

#' @export
bgnorm_results.QPTIFFImage <- function(x, ...) attr(x, "bgnorm_results")

#' @export
dim.QPTIFFImage <- function(x) {
  if (.is_lazy_qptiff(x)) dim(x$.da) else NextMethod()
}

#' @export
dimnames.QPTIFFImage <- function(x) {
  if (.is_lazy_qptiff(x)) dimnames(x$.da) else NextMethod()
}

#' Channel names of a QPTIFFImage (third-dimension names of the array)
#' @param x A \code{QPTIFFImage}.
#' @return Character vector of channel names, or \code{NULL}.
#' @export
names.QPTIFFImage <- function(x) dimnames(x)[[3L]]

#' Replace the dimnames of a QPTIFFImage
#'
#' Low-level setter used by \code{names<-.QPTIFFImage} and direct
#' \code{dimnames(img) <- } assignments.  Preserves the \code{QPTIFFImage}
#' class and the \code{"metadata"} attribute.
#'
#' @param x     A \code{QPTIFFImage}.
#' @param value A list of length 3 (\code{NULL} elements are allowed for the
#'   row and column dimensions).
#' @return \code{x} with updated dimnames.
#' @export
`dimnames<-.QPTIFFImage` <- function(x, value) {
  md <- attr(x, "metadata")
  br <- attr(x, "bgnorm_results")
  if (.is_lazy_qptiff(x)) {
    dimnames(x$.da) <- value
  } else {
    class(x) <- "array"
    dimnames(x) <- value
    class(x) <- c("QPTIFFImage", "array")
  }
  attr(x, "metadata") <- md
  attr(x, "bgnorm_results") <- br
  x
}

#' Set the channel names of a QPTIFFImage
#'
#' Replaces the third-dimension names (channel names) of a
#' \code{QPTIFFImage}.  Works for both eager (in-memory) and lazy (on-disk)
#' objects.
#'
#' @param x     A \code{QPTIFFImage}.
#' @param value Character vector of length equal to the number of channels.
#' @return \code{x} with updated channel names.
#' @export
#' @examples
#' arr <- array(1:8, dim = c(2, 2, 2),
#'              dimnames = list(NULL, NULL, c("ch1", "ch2")))
#' img <- as.QPTIFFImage(arr)
#' names(img) <- c("DAPI", "CD3")
#' names(img)
`names<-.QPTIFFImage` <- function(x, value) {
  nch <- dim(x)[3L]
  if (length(value) != nch)
    stop("Length of replacement (", length(value), ") must equal the number ",
         "of channels (", nch, ").")
  dn       <- dimnames(x)
  if (is.null(dn)) dn <- vector("list", 3L)
  dn[[3L]] <- as.character(value)
  dimnames(x) <- dn
  x
}

#' Number of channels in a QPTIFFImage
#' @param x A \code{QPTIFFImage}.
#' @return Integer; the number of channels.
#' @export
length.QPTIFFImage <- function(x) dim(x)[3L]

#' Subset a QPTIFFImage
#'
#' Standard 3-D array subscripting with all missing-subscript combinations
#' supported.  The result is \emph{always} a \code{QPTIFFImage}, including
#' when a single channel is selected.  To extract a plain 2-D matrix for one
#' channel, use \code{unclass(img[, , "DAPI"])[, , 1L]}.
#' \itemize{
#'   \item \code{img[, , "DAPI"]} - single channel (returns a 1-channel QPTIFFImage)
#'   \item \code{img[, , c("DAPI","CD3")]} - multi-channel sub-image
#'   \item \code{img[1:512, 1:512, ]} - spatial crop
#' }
#' @param x    A \code{QPTIFFImage}.
#' @param i    Row indices (spatial), or missing.
#' @param j    Column indices (spatial), or missing.
#' @param k    Channel indices (numeric or character), or missing.
#' @param drop Ignored; the result is always a 3-D \code{QPTIFFImage}.
#' @return A \code{QPTIFFImage}.
#' @export
`[.QPTIFFImage` <- function(x, i, j, k, drop = FALSE) {
  md           <- attr(x, "metadata")
  br           <- attr(x, "bgnorm_results")
  # bgnorm_results are per-pixel vectors - they are only valid when the spatial
  # dimensions are untouched.  Spatial subsetting (i or j specified) invalidates
  # them; channel-only subsetting (only k specified, or nothing) preserves them.
  spatial_sub  <- !missing(i) || !missing(j)

  .subset_br <- function(arr3d) {
    if (is.null(br) || spatial_sub) return(NULL)
    sel <- dimnames(arr3d)[[3L]]
    out <- br[intersect(sel, names(br))]
    if (length(out) == 0L) NULL else out
  }

  if (.is_lazy_qptiff(x)) {
    da <- x$.da
    res <- if (missing(i) && missing(j) && missing(k)) {
      da[, , , drop = FALSE]
    } else if (missing(i) && missing(j)) {
      da[, , k, drop = FALSE]
    } else if (missing(i) && missing(k)) {
      da[, j, , drop = FALSE]
    } else if (missing(j) && missing(k)) {
      da[i, , , drop = FALSE]
    } else if (missing(i)) {
      da[, j, k, drop = FALSE]
    } else if (missing(j)) {
      da[i, , k, drop = FALSE]
    } else if (missing(k)) {
      da[i, j, , drop = FALSE]
    } else {
      da[i, j, k, drop = FALSE]
    }
    arr <- as.array(res)
    if (length(dim(arr)) == 3L)
      .new_QPTIFFImage(arr, md, bgnorm_results = .subset_br(arr))
    else arr
  } else {
    class(x) <- "array"
    # Always use drop = FALSE so the result is 3-D and can be wrapped in
    # QPTIFFImage.  We handle every combination of missing subscripts explicitly
    # because missing symbols do not propagate reliably across all R versions
    # when passed positionally to a primitive like `[`.
    res <- if (missing(i) && missing(j) && missing(k)) {
      x[, , , drop = FALSE]
    } else if (missing(i) && missing(j)) {
      x[, , k, drop = FALSE]
    } else if (missing(i) && missing(k)) {
      x[, j, , drop = FALSE]
    } else if (missing(j) && missing(k)) {
      x[i, , , drop = FALSE]
    } else if (missing(i)) {
      x[, j, k, drop = FALSE]
    } else if (missing(j)) {
      x[i, , k, drop = FALSE]
    } else if (missing(k)) {
      x[i, j, , drop = FALSE]
    } else {
      x[i, j, k, drop = FALSE]
    }
    if (is.array(res) && length(dim(res)) == 3L)
      .new_QPTIFFImage(res, md, bgnorm_results = .subset_br(res))
    else
      res
  }
}

#' Convert a QPTIFFImage to a named list of 2-D channel matrices
#' @param x   A \code{QPTIFFImage}.
#' @param ... Unused.
#' @return A named list of 2-D matrices, one per channel.
#' @export
as.list.QPTIFFImage <- function(x, ...) {
  chs <- names(x)
  n   <- length(chs)
  if (.is_lazy_qptiff(x)) {
    da  <- x$.da
    out <- lapply(seq_len(n), function(k) as.array(da[, , k, drop = FALSE])[, , 1L])
  } else {
    arr <- unclass(x)
    out <- lapply(seq_len(n), function(k) arr[, , k])
  }
  names(out) <- chs
  out
}

#' Coerce a 3-D array to a QPTIFFImage
#'
#' Attaches the \code{QPTIFFImage} class to any 3-D numeric array
#' \code{[H x W x C]}.  Channel names are taken from \code{dimnames(x)[[3]]}.
#' This is the intended way to promote a plain array produced outside
#' \code{read_qptiff} (e.g., the adjusted output of \code{bgnorm_markers})
#' into a first-class \code{QPTIFFImage}.
#'
#' @param x        A 3-D numeric array \code{[H x W x C]}.  For a
#'   \code{QPTIFFImage}, the object is returned unchanged.
#' @param metadata Optional named list of metadata to embed.  Ignored when
#'   \code{x} is already a \code{QPTIFFImage}.
#' @param ...      Unused; present for S3 method consistency.
#'
#' @return A \code{QPTIFFImage}.
#' @export
#' @examples
#' arr <- array(runif(20 * 20 * 3), dim = c(20, 20, 3),
#'              dimnames = list(NULL, NULL, c("DAPI", "CD3", "CD8")))
#' img <- as.QPTIFFImage(arr)
#' class(img)
as.QPTIFFImage <- function(x, ...) UseMethod("as.QPTIFFImage")

#' @rdname as.QPTIFFImage
#' @export
as.QPTIFFImage.QPTIFFImage <- function(x, ...) x

#' @rdname as.QPTIFFImage
#' @export
as.QPTIFFImage.array <- function(x, metadata = list(), ...) {
  if (length(dim(x)) != 3L)
    stop("'x' must be a 3-D array [H x W x C]; got dim: ",
         paste(dim(x), collapse = " x "))
  .new_QPTIFFImage(x, metadata = metadata)
}

#' @rdname as.QPTIFFImage
#' @export
as.QPTIFFImage.DelayedArray <- function(x, metadata = list(), ...) {
  .new_QPTIFFImage_lazy(x, metadata = metadata)
}

#' Materialise a QPTIFFImage to a plain 3-D array
#'
#' For eager QPTIFFImages, equivalent to \code{unclass(x)}.  For lazy
#' (on-disk) QPTIFFImages, reads all pixel data from disk.
#'
#' @param x   A \code{QPTIFFImage}.
#' @param ... Unused.
#' @return A plain numeric array of dimension \code{c(H, W, C)}.
#' @export
as.array.QPTIFFImage <- function(x, ...) {
  if (.is_lazy_qptiff(x)) as.array(x$.da) else unclass(x)
}

#' Coerce a single-channel QPTIFFImage to a matrix
#'
#' Drops the singleton channel dimension and returns a plain \code{H x W}
#' matrix.  The channel name is stored in \code{attr(result, "channel")}.
#' If the object has more than one channel, use \code{img[, , "DAPI"]} to
#' select one first.
#'
#' @param x   A \code{QPTIFFImage} with exactly one channel, \emph{or} a
#'   2-D array (already a matrix-like object), which is returned with
#'   \code{\link[base]{as.matrix}} applied directly.
#' @param ... Unused.
#' @return A numeric matrix of dimension \code{c(H, W)}.
#' @export
#' @examples
#' arr <- array(1:12, dim = c(3, 4, 1),
#'              dimnames = list(NULL, NULL, "DAPI"))
#' img <- as.QPTIFFImage(arr)
#' m   <- as.matrix(img)
#' dim(m)          # 3 x 4
#' attr(m, "channel")  # "DAPI"
as.matrix.QPTIFFImage <- function(x, ...) {
  d <- dim(x)
  if (length(d) == 2L)
    return(NextMethod())
  if (d[3L] != 1L)
    stop("as.matrix() requires a single-channel QPTIFFImage (got ", d[3L],
         " channels).\nSubset first: img[, , \"", names(x)[1L], "\"]")
  ch  <- names(x)[1L]
  mat <- if (.is_lazy_qptiff(x))
    drop(as.array(x$.da[, , 1L, drop = FALSE]))
  else
    unclass(x)[, , 1L]
  if (!is.null(ch)) attr(mat, "channel") <- ch
  mat
}

#' Access the metadata embedded in a QPTIFFImage
#'
#' Retrieves the \code{\link{QPTIFFMetadata}} object (an OME-organised
#' hierarchy of slide / image / channel / scale metadata) stored inside a
#' \code{\link{QPTIFFImage}}.  The generic falls back to
#' \code{\link[S4Vectors]{metadata}} for S4 objects such as
#' \code{SummarizedExperiment}, so \code{metadata(img)} and
#' \code{metadata(spe)} both work when bgnormR is loaded.
#'
#' The returned object mirrors the OME hierarchy:
#' \describe{
#'   \item{\code{slide}}{Slide-constant identity (id, barcode, operator, ...).}
#'   \item{\code{images}}{One scene per distinct image, each with
#'     \code{image_info} (optics, camera, scan_resolution), \code{channels}
#'     (per-channel \code{name}, \code{fluorophore}, \code{dye_from_name},
#'     \code{exposure_time_us}, wavelengths, \code{color_rgb}, ...) and
#'     \code{scales} (per pyramid level).}
#'   \item{\code{acquisition_format}}{One of \code{"brightfield"},
#'     \code{"polaris_scanband"}, \code{"fusion_paged"}, \code{"ome_tiff"},
#'     \code{"ome_zarr"}, or \code{"unknown"}.}
#' }
#' Use the accessors (\code{\link{qpi_format}}, \code{\link{qpi_channels}},
#' \code{\link{qpi_pixel_size_um}}, \code{\link{channel_table}}, ...) rather
#' than indexing the nested list directly.
#'
#' @param x A \code{\link{QPTIFFImage}}.
#' @param ... Unused.
#' @return A \code{\link{QPTIFFMetadata}} object (a named list).
#'
#' @export
#' @examples
#' path <- system.file("extdata", "PA_HNC_sample.ome.tiff", package = "bgnormR")
#' img  <- read_qptiff(path)
#' meta <- metadata(img)
#' qpi_format(meta)
#' # Per-channel dye information as a tidy table
#' channel_table(meta)[, c("name", "fluorophore", "dye_from_name")]
metadata <- function(x, ...) UseMethod("metadata")

#' @rdname metadata
#' @export
metadata.QPTIFFImage <- function(x, ...) attr(x, "metadata")

# Fallback for S4 objects (SummarizedExperiment, etc.) - delegates to S4Vectors.
#' @export
metadata.default <- function(x, ...) {
  if (isS4(x))
    S4Vectors::metadata(x, ...)
  else
    stop("No 'metadata' method for class '", class(x)[1L], "'.")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ============================================================
# Protein / dye name splitting
# ============================================================

# End-anchored, case-insensitive regex matching common fluorescence dye suffixes
# that are appended to protein names in compound biomarker tags (e.g. "CD3e-AF647").
.DYE_SUFFIX_RE <- paste0(
  "(?i)-(",
  "AF[0-9]+",                       # Alexa Fluor: AF647, AF750, AF488, AF405
  "|Atto[0-9]+[NW]?",               # Atto dyes: Atto550, Atto647N
  "|OPAL[0-9]+",                    # Polaris OPAL: OPAL520, OPAL690
  "|Cy[0-9]+(?:[.][0-9]+)?",        # Cyanine: Cy3, Cy5, Cy5.5
  "|PE(?:-Cy[0-9]+)?",              # PE, PE-Cy5, PE-Cy7
  "|APC(?:-Cy[0-9]+)?",             # APC, APC-Cy7
  "|FITC",                          # FITC
  "|BV[0-9]+",                      # Brilliant Violet: BV421, BV786
  "|PerCP(?:-Cy[0-9.]+)?",          # PerCP, PerCP-Cy5.5
  "|BUV[0-9]+",                     # Brilliant UV: BUV395, BUV737
  "|Pacific(?:Blue|Orange)",        # Pacific Blue / Pacific Orange
  ")$"
)

#' Split a compound biomarker name into protein and fluorescence dye
#'
#' Parses names of the form \code{"CD3e-AF647"} into a protein component
#' (\code{"CD3e"}) and a fluorescence dye component (\code{"AF647"}).
#' The already-parsed \code{fluorophore} field is used as the primary signal;
#' a regex against known dye naming conventions serves as a fallback.
#'
#' Handles multi-hyphen names correctly: \code{"HLA-A-Atto550"} ->
#' protein \code{"HLA-A"}, dye \code{"Atto550"}.
#'
#' @param name       Character scalar; raw channel / biomarker name.
#' @param fluorophore Character scalar or \code{NULL}; fluorophore already
#'   parsed from the XML \code{<Fluorophore>} element.
#' @return A list with elements \code{protein} (character) and \code{dye}
#'   (character or \code{NULL}).
#' @keywords internal
.split_protein_dye <- function(name, fluorophore = NULL) {
  if (!is.character(name) || !nzchar(name))
    return(list(protein = name, dye = NULL))

  # Primary: if the parsed fluorophore is known and the name ends with
  # -{fluorophore} (case-insensitive), use that to locate the split point.
  # Return the dye with the casing as it appears in the name, not the tag.
  if (!is.null(fluorophore) && nzchar(fluorophore)) {
    suffix <- paste0("-", fluorophore)
    slen   <- nchar(suffix)
    nlen   <- nchar(name)
    if (nlen > slen &&
        toupper(substr(name, nlen - slen + 1L, nlen)) == toupper(suffix)) {
      return(list(
        protein = substr(name, 1L, nlen - slen),
        dye     = substr(name, nlen - slen + 2L, nlen)   # +2 to skip the "-"
      ))
    }
  }

  # Fallback: match end-anchored dye pattern in the name string.
  m <- regmatches(name, regexpr(.DYE_SUFFIX_RE, name, perl = TRUE))
  if (length(m) == 1L && nzchar(m))
    return(list(
      protein = substr(name, 1L, nchar(name) - nchar(m)),
      dye     = sub("^-", "", m)
    ))

  list(protein = name, dye = NULL)
}

# ============================================================
# Main QPTIFF XML parser  (ported from qptiff_parser.py)
# ============================================================

#' Parse a PerkinElmer QPI ImageDescription XML
#'
#' @param xml_str       Raw XML string from TIFF tag 270 of the first page.
#' @param per_page_xmls Character vector of XML strings, one per channel page.
#' @param n_channels    Fallback channel count when XML is absent.
#' @param is_rgb        TRUE when the first page is RGB (brightfield signal).
#'
#' @return Named list with \code{slide}, \code{image_info}, \code{channels},
#'   \code{format}, \code{raw_xml}.
#' @keywords internal
.parse_qpi_xml <- function(xml_str, per_page_xmls = NULL,
                            n_channels = NULL, is_rgb = FALSE) {
  slide      <- list()
  image_info <- list()
  channels   <- list()

  empty <- function() list(slide = .new_slide_info(), image_info = .new_image_info(),
                           channels = channels, format = .QPTIFF_UNK,
                           raw_xml = xml_str %||% "")

  if (is.null(xml_str) || !nzchar(xml_str)) return(empty())

  doc <- tryCatch(xml2::read_xml(xml_str), error = function(e) NULL)
  if (is.null(doc)) {
    message("Failed to parse QPI XML.")
    return(empty())
  }

  root_name <- xml2::xml_name(doc)
  if (!grepl("PerkinElmer.QPI|PerkinElmerQPI", root_name, perl = TRUE))
    message("Note: XML root '", root_name, "' not a recognised QPI description.")

  # ---- Slide-level fields -----------------------------------------------
  slide <- list(
    description_version  = .xt(doc, "DescriptionVersion"),
    acquisition_software = .xt(doc, "AcquisitionSoftware"),
    identifier           = .xt(doc, "Identifier"),
    slide_id             = .xt(doc, "SlideID"),
    barcode              = .xt(doc, "Barcode"),
    study_name           = .xt(doc, "StudyName"),
    operator_name        = .xt(doc, "OperatorName"),
    computer_name        = .xt(doc, "ComputerName"),
    instrument_type      = .xt(doc, "InstrumentType"),
    validation_code      = .xt(doc, "ValidationCode"),
    sample_description   = .xt(doc, "SampleDescription")
  )

  # ---- Image-level fields -----------------------------------------------
  image_info <- list(
    image_type  = .xt(doc, "ImageType"),
    bf_lamp_type = .xt(doc, "BFLampType"),
    lamp_type   = .xt(doc, "LampType"),
    objective   = .xt(doc, "Objective")
  )

  # ScanProfile
  sp <- xml2::xml_find_first(doc, "ScanProfile")
  if (!inherits(sp, "xml_missing")) {
    # Polaris format nests an extra <root> inside <ScanProfile>
    sp_root_candidate <- xml2::xml_find_first(sp, "root")
    sp_root <- if (inherits(sp_root_candidate, "xml_missing")) sp else sp_root_candidate

    image_info$scan_profile_name <- .xt(sp_root, "Name")
    image_info$scan_mode         <- .xt(sp_root, "Mode")
    image_info$is_tma            <- .xb(sp_root, "SampleIsTMA")
    image_info$opal_kit_type     <- .xt(sp_root, "OpalKitType")
    image_info$compression       <- .xt(sp_root, "Compression")
    image_info$jpeg_quality      <- .xi(sp_root, "JPEGQuality")
    image_info$saturation_protection_type <- .xt(sp_root, "SaturationProtectionType") %||%
                                             .xt(sp_root, "FieldSaturationProtectionType")
    image_info$coverslip_thickness <- .xt(sp_root, "CoverslipThickness")
    image_info$is_rna            <- .xb(sp_root, "IsRNA")
    cam_sp <- xml2::xml_find_first(sp_root, "CameraSettings")
    if (!inherits(cam_sp, "xml_missing")) {
      image_info$rotate_image  <- .xb(cam_sp, "RotateImage")
      image_info$mirror_image  <- .xb(cam_sp, "MirrorImage")
    }

    # Fusion 1.x: ScanProfile is a JSON blob
    sp_text <- trimws(xml2::xml_text(sp))
    if (startsWith(sp_text, "{")) {
      image_info <- .parse_scan_profile_json(image_info, sp_text)
    }
  }

  image_info <- c(image_info, .parse_scan_resolution_xml(doc))
  base_px <- image_info$pixel_size_um
  if (!is.null(base_px)) {
    image_info$scale_factor      <- base_px
    image_info$scale_factor_unit <- "um"
  }
  image_info$camera <- .parse_camera_xml(doc)

  exposure_times <- .parse_exposure_times_xml(doc)
  scan_mode      <- image_info$scan_mode %||% ""

  # ---- Format detection -------------------------------------------------
  fmt <- .detect_qptiff_format(doc, scan_mode, is_rgb)

  # ---- Channel parsing per format ---------------------------------------
  if (fmt == .QPTIFF_BF) {
    n <- if (!is.null(n_channels) && n_channels > 0L) n_channels else 3L
    channels <- .parse_brightfield_channels(n, doc)

  } else if (fmt == .QPTIFF_PSB) {
    distinct <- length(unique(Filter(nzchar, per_page_xmls %||% character(0))))
    if (distinct > 1L && !is.null(per_page_xmls)) {
      channels <- .parse_per_page_channels(per_page_xmls)
    } else {
      channels <- .parse_polaris_channels(doc, exposure_times)
    }

  } else if (fmt == .QPTIFF_FP) {
    if (!is.null(per_page_xmls) && length(per_page_xmls) > 0L) {
      channels <- .parse_per_page_channels(per_page_xmls)
    } else if (!is.null(n_channels) && n_channels > 0L) {
      channels <- lapply(seq_len(n_channels) - 1L, function(i) {
        et <- if (i < length(exposure_times)) exposure_times[[i + 1L]] else NULL
        list(index = i, name = paste0("Channel_", i),
             exposure_time_us = et, is_brightfield = FALSE)
      })
    }

  } else {  # UNKNOWN: try strategies in order
    flu <- .parse_polaris_channels(doc, exposure_times)
    if (length(flu) > 0L) {
      channels <- flu
    } else if (!is.null(per_page_xmls) && length(per_page_xmls) > 0L) {
      channels <- .parse_per_page_channels(per_page_xmls)
    } else if (!is.null(n_channels) && n_channels > 0L) {
      channels <- lapply(seq_len(n_channels) - 1L, function(i) {
        et <- if (i < length(exposure_times)) exposure_times[[i + 1L]] else NULL
        list(index = i, name = paste0("Channel_", i),
             exposure_time_us = et, is_brightfield = FALSE)
      })
    }
  }

  # Fill in missing exposure times
  for (i in seq_along(channels)) {
    if (is.null(channels[[i]]$exposure_time_us)) {
      idx <- channels[[i]]$index + 1L
      if (idx <= length(exposure_times))
        channels[[i]]$exposure_time_us <- exposure_times[[idx]]
    }
  }

  # Fusion paged: fluorophore from JSON ScanProfile
  if (fmt == .QPTIFF_FP && any(vapply(channels, function(ch) is.null(ch$fluorophore), logical(1L)))) {
    sp_elem <- xml2::xml_find_first(doc, "ScanProfile")
    if (!inherits(sp_elem, "xml_missing")) {
      sp_txt <- trimws(xml2::xml_text(sp_elem))
      if (startsWith(sp_txt, "{"))
        channels <- .supplement_fluorophore_from_json(channels, sp_txt)
    }
  }

  list(
    slide      = do.call(.new_slide_info, slide),
    image_info = .nest_image_info(image_info),
    channels   = channels,
    format     = fmt,
    raw_xml    = xml_str
  )
}

# ---- Reorganise a flat QPI image_info into the nested ImageInfo schema -----
#
# The QPI parser accumulates image-level fields in a flat list; here we route
# them into the canonical ImageInfo { camera, scan_resolution } shape.  Fields
# with no schema home (lamp_type, compression, coverslip_thickness, ...) are
# preserved under image_info$extra so nothing is silently dropped.

#' @keywords internal
.nest_image_info <- function(flat) {
  cam <- flat$camera %||% list()
  camera <- .new_camera_info(
    camera_name = cam$camera_name, camera_type = cam$camera_type,
    gain = cam$gain, bit_depth = cam$bit_depth, binning = cam$binning
  )
  scan_res <- .new_scan_resolution(
    magnification      = flat$magnification,
    objective_name     = flat$objective_name,
    binning            = flat$binning,
    base_pixel_size_um = flat$pixel_size_um
  )
  ii <- .new_image_info(
    image_type        = flat$image_type,
    objective         = flat$objective,
    bf_lamp_type      = flat$bf_lamp_type,
    scan_profile_name = flat$scan_profile_name,
    scan_mode         = flat$scan_mode,
    is_tma            = flat$is_tma,
    opal_kit_type     = flat$opal_kit_type,
    xposition_um      = flat$xposition_um,
    yposition_um      = flat$yposition_um,
    camera            = camera,
    scan_resolution   = scan_res
  )
  used   <- c("image_type", "objective", "bf_lamp_type", "scan_profile_name",
              "scan_mode", "is_tma", "opal_kit_type", "xposition_um",
              "yposition_um", "camera", "magnification", "objective_name",
              "binning", "pixel_size_um", "scale_factor", "scale_factor_unit")
  extras <- flat[setdiff(names(flat), used)]
  extras <- extras[!vapply(extras, is.null, logical(1L))]
  if (length(extras) > 0L) ii$extra <- extras
  ii
}

# ---- Format detection  (== _detect_format) --------------------------------

#' @keywords internal
.detect_qptiff_format <- function(doc, scan_mode, is_rgb) {
  if (grepl("brightfield", scan_mode, ignore.case = TRUE) || is_rgb)
    return(.QPTIFF_BF)

  if (length(xml2::xml_find_all(doc, ".//ScanBands-i")) > 0L)
    return(.QPTIFF_PSB)

  sp <- xml2::xml_find_first(doc, "ScanProfile")
  if (!inherits(sp, "xml_missing")) {
    sp_text <- trimws(xml2::xml_text(sp))
    if (startsWith(sp_text, "{"))
      return(.QPTIFF_FP)
  }

  .QPTIFF_UNK
}

# ---- Exposure time parsing (== _parse_exposure_times) --------------------

#' @keywords internal
.parse_exposure_times_xml <- function(doc) {
  eta <- xml2::xml_find_first(doc, "ExposureTimeArray")
  if (!inherits(eta, "xml_missing")) {
    vals <- xml2::xml_find_all(eta, "Value")
    return(lapply(vals, function(v) .xf_node(v)))
  }
  et <- xml2::xml_find_first(doc, "ExposureTime")
  list(if (inherits(et, "xml_missing")) NULL else .xf_node(et))
}

# ---- Camera parsing (== _parse_camera) -----------------------------------

#' @keywords internal
.parse_camera_xml <- function(doc) {
  cs <- xml2::xml_find_first(doc, "CameraSettings")
  list(
    camera_name = .xt(doc, "CameraName"),
    camera_type = .xt(doc, "CameraType"),
    gain        = if (inherits(cs, "xml_missing")) NULL else .xf(cs, "Gain"),
    bit_depth   = if (inherits(cs, "xml_missing")) NULL else .xi(cs, "BitDepth"),
    binning     = if (inherits(cs, "xml_missing")) NULL else .xi(cs, "Binning")
  )
}

# ---- Scan resolution (== _parse_scan_resolution) -------------------------

#' @keywords internal
.parse_scan_resolution_xml <- function(doc) {
  sr <- xml2::xml_find_first(doc, ".//ScanResolution")
  if (inherits(sr, "xml_missing"))
    return(list(pixel_size_um = NULL, magnification = NULL,
                objective_name = NULL, binning = NULL))
  list(
    pixel_size_um  = .xf(sr, "PixelSizeMicrons"),
    magnification  = .xf(sr, "Magnification"),
    objective_name = .xt(sr, "ObjectiveName"),
    binning        = .xi(sr, "Binning")
  )
}

# ---- Active band selection (== _select_active_band) ----------------------

#' @keywords internal
.select_active_band <- function(filter_node) {
  if (inherits(filter_node, "xml_missing") || is.null(filter_node)) return(NULL)
  bands <- xml2::xml_find_all(filter_node, "Bands/Band")
  if (length(bands) == 0L) return(NULL)
  for (b in bands) {
    act <- xml2::xml_find_first(b, "Active")
    if (!inherits(act, "xml_missing")) {
      txt <- tolower(trimws(xml2::xml_text(act)))
      if (txt %in% c("true", "1", "yes")) return(b)
    }
  }
  bands[[1L]]
}

# ---- Per-channel field extraction (== _populate_channel_fields_from_element) --

#' @keywords internal
.populate_channel_from_node <- function(node, fields) {
  # Biomarker name (priority order from Python _BIOMARKER_TAGS)
  if (is.null(fields$name) || startsWith(fields$name, "Channel_")) {
    for (tag in .BIOMARKER_TAGS) {
      el <- xml2::xml_find_first(node, tag)
      if (!inherits(el, "xml_missing")) {
        val <- trimws(xml2::xml_text(el))
        if (nzchar(val) && !val %in% c("None", "--")) {
          fields$name <- val
          break
        }
      }
    }
  }

  # Fluorophore
  if (is.null(fields$fluorophore)) {
    for (tag in .FLUOROPHORE_TAGS) {
      el <- xml2::xml_find_first(node, tag)
      if (!inherits(el, "xml_missing")) {
        val <- trimws(xml2::xml_text(el))
        if (nzchar(val)) { fields$fluorophore <- val; break }
      }
    }
  }

  # Exposure time (us)
  et_node <- xml2::xml_find_first(node, "ExposureTime")
  if (!inherits(et_node, "xml_missing") && is.null(fields$exposure_time_us))
    fields$exposure_time_us <- .xf_node(et_node)

  # Emission wavelength via active band
  em_filt  <- xml2::xml_find_first(node, ".//EmissionFilter")
  em_band  <- .select_active_band(em_filt)
  if (!is.null(em_band) && is.null(fields$emission_wavelength_nm)) {
    cuton  <- .xf(em_band, "Cuton")
    cutoff <- .xf(em_band, "Cutoff")
    if (!is.null(cuton) && !is.null(cutoff))
      fields$emission_wavelength_nm <- (cuton + cutoff) / 2
  }

  # Excitation wavelength via active band
  ex_filt <- xml2::xml_find_first(node, ".//ExcitationFilter")
  ex_band <- .select_active_band(ex_filt)
  if (!is.null(ex_band) && is.null(fields$excitation_wavelength_nm)) {
    cuton  <- .xf(ex_band, "Cuton")
    cutoff <- .xf(ex_band, "Cutoff")
    if (!is.null(cuton) && !is.null(cutoff))
      fields$excitation_wavelength_nm <- (cuton + cutoff) / 2
  }

  # Fallback wavelengths (old ScanBand format)
  if (is.null(fields$emission_wavelength_nm))
    fields$emission_wavelength_nm  <- .xf(node, ".//HomeWavelength")
  if (is.null(fields$excitation_wavelength_nm))
    fields$excitation_wavelength_nm <- .xf(node, ".//ExcitationWavelength")

  # Infer emission from OPAL fluorophore name (e.g. "OPAL520" -> 520)
  if (is.null(fields$emission_wavelength_nm) && !is.null(fields$fluorophore)) {
    m <- regmatches(fields$fluorophore, regexpr("[0-9]{3,4}", fields$fluorophore))
    if (length(m) == 1L) fields$emission_wavelength_nm <- as.numeric(m)
  }

  # Display colour <Color>r,g,b</Color>
  if (is.null(fields$color_rgb)) {
    col_node <- xml2::xml_find_first(node, "Color")
    if (inherits(col_node, "xml_missing"))
      col_node <- xml2::xml_find_first(node, ".//Color")
    fields$color_rgb <- .parse_color_node(col_node)
  }

  fields$is_unmixed_component <- .xb(node, "IsUnmixedComponent") %||%
                                  fields$is_unmixed_component
  fields$signal_units         <- .xi(node, "SignalUnits") %||% fields$signal_units
  fields$objective            <- .xt(node, "Objective") %||% fields$objective
  fields$autofluorescence_subtracted <- .xb(node, "AutofluorescenceSubtracted") %||%
                                         fields$autofluorescence_subtracted
  fields$auto_expose_type     <- .xt(node, ".//AutoExposeType") %||% fields$auto_expose_type

  # Responsivity
  resp <- xml2::xml_find_first(node, "Responsivity")
  if (!inherits(resp, "xml_missing")) {
    rf <- xml2::xml_find_first(resp, "Filter")
    if (!inherits(rf, "xml_missing")) {
      fields$responsivity               <- .xf(rf, "Response")
      fields$responsivity_filter_id     <- .xt(rf, "FilterID")
      fields$responsivity_date          <- .xt(rf, "Date")
      fields$responsivity_filter_name   <- .xt(rf, "Name")
    }
  }

  # Filter identification
  ex_f <- xml2::xml_find_first(node, "ExcitationFilter")
  if (!inherits(ex_f, "xml_missing")) {
    fields$excitation_filter_name         <- .xt(ex_f, "Name")
    fields$excitation_filter_manufacturer <- .xt(ex_f, "Manufacturer")
    fields$excitation_filter_part_no      <- .xt(ex_f, "PartNo")
  }
  em_f <- xml2::xml_find_first(node, "EmissionFilter")
  if (!inherits(em_f, "xml_missing")) {
    fields$emission_filter_name         <- .xt(em_f, "Name")
    fields$emission_filter_manufacturer <- .xt(em_f, "Manufacturer")
    fields$emission_filter_part_no      <- .xt(em_f, "PartNo")
  }

  # Per-channel camera settings
  cs <- xml2::xml_find_first(node, "CameraSettings")
  if (!inherits(cs, "xml_missing")) {
    fields$gain             <- .xf(cs, "Gain") %||% fields$gain
    fields$binning          <- .xi(cs, "Binning") %||% fields$binning
    fields$bit_depth        <- .xi(cs, "BitDepth") %||% fields$bit_depth
    fields$offset_counts    <- .xi(cs, "OffsetCounts") %||% fields$offset_counts
    fields$camera_orientation <- .xt(cs, "Orientation") %||% fields$camera_orientation
    roi <- xml2::xml_find_first(cs, "ROI")
    if (!inherits(roi, "xml_missing")) {
      fields$roi_x      <- .xi(roi, "X")
      fields$roi_y      <- .xi(roi, "Y")
      fields$roi_width  <- .xi(roi, "Width")
      fields$roi_height <- .xi(roi, "Height")
    }
  } else {
    # Older ScanBands: Gain/Binning may be direct children
    fields$gain    <- .xf(node, ".//Gain") %||% fields$gain
    fields$binning <- .xi(node, ".//Binning") %||% fields$binning
  }

  fields
}

# ---- Brightfield channel parsing (== _parse_brightfield_channels) ---------

#' @keywords internal
.parse_brightfield_channels <- function(n_samples, doc = NULL) {
  nms <- if (n_samples == 3L) c("Red", "Green", "Blue") else
          paste0("Sample_", seq_len(n_samples) - 1L)
  shared <- list()
  if (!is.null(doc)) {
    et_node <- xml2::xml_find_first(doc, "ExposureTime")
    if (!inherits(et_node, "xml_missing"))
      shared$exposure_time_us <- .xf_node(et_node)
    shared$signal_units         <- .xi(doc, "SignalUnits")
    shared$is_unmixed_component <- .xb(doc, "IsUnmixedComponent")
    shared$objective            <- .xt(doc, "Objective")
    cs <- xml2::xml_find_first(doc, "CameraSettings")
    if (!inherits(cs, "xml_missing")) {
      shared$gain    <- .xf(cs, "Gain")
      shared$binning <- .xi(cs, "Binning")
      shared$bit_depth <- .xi(cs, "BitDepth")
      shared$offset_counts <- .xi(cs, "OffsetCounts")
      shared$camera_orientation <- .xt(cs, "Orientation")
      roi <- xml2::xml_find_first(cs, "ROI")
      if (!inherits(roi, "xml_missing")) {
        shared$roi_x      <- .xi(roi, "X")
        shared$roi_y      <- .xi(roi, "Y")
        shared$roi_width  <- .xi(roi, "Width")
        shared$roi_height <- .xi(roi, "Height")
      }
    }
  }
  lapply(seq_len(n_samples), function(i) {
    c(list(index = i - 1L, name = nms[i], is_brightfield = TRUE), shared)
  })
}

# ---- Polaris ScanBand channel parsing (== _parse_fluorescence_channels) ---

#' @keywords internal
.parse_polaris_channels <- function(doc, exposure_times) {
  bands <- xml2::xml_find_all(doc, ".//ScanBands-i")
  if (length(bands) == 0L) return(list())
  lapply(seq_along(bands), function(i) {
    band   <- bands[[i]]
    fields <- list(index = i - 1L, name = paste0("Channel_", i - 1L),
                   is_brightfield = FALSE)
    fields <- .populate_channel_from_node(band, fields)
    if (is.null(fields$exposure_time_us) && i <= length(exposure_times))
      fields$exposure_time_us <- exposure_times[[i]]
    fields
  })
}

# ---- Fusion per-page channel parsing (== _parse_channels_from_per_page_xmls) --

#' @keywords internal
.parse_per_page_channels <- function(per_page_xmls) {
  lapply(seq_along(per_page_xmls), function(i) {
    xml_str <- per_page_xmls[[i]]
    fields  <- list(index = i - 1L, name = paste0("Channel_", i - 1L),
                    is_brightfield = FALSE)
    if (!is.null(xml_str) && nzchar(xml_str)) {
      page_doc <- tryCatch(xml2::read_xml(xml_str), error = function(e) NULL)
      if (!is.null(page_doc))
        fields <- .populate_channel_from_node(page_doc, fields)
    }
    fields
  })
}

# ---- Fluorophore supplement from Fusion JSON (== inline logic in parse_qpi_xml) --

#' @keywords internal
.supplement_fluorophore_from_json <- function(channels, json_str) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(channels)
  sp <- tryCatch(jsonlite::fromJSON(json_str, simplifyVector = FALSE),
                 error = function(e) NULL)
  if (is.null(sp) || !is.list(sp)) return(channels)
  filter_fluors <- tryCatch({
    ch_list <- sp$experimentDescription$channels
    if (is.list(ch_list)) vapply(ch_list, function(c) c$name %||% "", character(1L))
    else character(0L)
  }, error = function(e) character(0L))
  filter_fluors <- Filter(nzchar, filter_fluors)
  if (length(filter_fluors) == 0L) return(channels)
  n <- length(filter_fluors)
  for (i in seq_along(channels)) {
    if (is.null(channels[[i]]$fluorophore))
      channels[[i]]$fluorophore <- filter_fluors[((channels[[i]]$index) %% n) + 1L]
  }
  channels
}

# ---- Scan profile JSON (== _parse_scan_profile_json) ----------------------

#' @keywords internal
.parse_scan_profile_json <- function(image_info, json_str) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(image_info)
  sp <- tryCatch(jsonlite::fromJSON(json_str, simplifyVector = FALSE),
                 error = function(e) NULL)
  if (is.null(sp) || !is.list(sp)) return(image_info)
  if (is.null(image_info$is_tma) && is.logical(sp$isTma))
    image_info$is_tma <- sp$isTma
  if (is.null(image_info$binning) && is.numeric(sp$binning))
    image_info$binning <- as.integer(sp$binning)
  exp <- sp$experimentDescription
  if (is.list(exp) && is.null(image_info$scan_profile_name)) {
    nm <- exp$name
    if (is.character(nm) && nzchar(nm)) image_info$scan_profile_name <- nm
  }
  image_info
}

# ============================================================
# IFD chain reader - reads ImageDescription + SamplesPerPixel
# from every TIFF/BigTIFF page
# ============================================================

#' Read tag 270 and tag 277 from every IFD in a TIFF/BigTIFF file
#'
#' @param path      File path.
#' @param max_pages Maximum IFDs to traverse (default 1000).
#'
#' @return List with \code{descriptions} (character) and
#'   \code{samples_per_pixel} (integer).
#' @keywords internal
.read_all_ifd_info <- function(path, max_pages = 1000L) {
  f <- tryCatch(file(path, "rb"), error = function(e) NULL)
  if (is.null(f)) return(list(descriptions = character(0L),
                               samples_per_pixel = integer(0L)))
  on.exit(close(f), add = TRUE)

  # Header
  bom    <- readBin(f, "raw", n = 2L)
  endian <- if (identical(bom, as.raw(c(0x49L, 0x49L)))) "little" else "big"
  magic  <- readBin(f, "integer", n = 1L, size = 2L, signed = FALSE,
                    endian = endian)
  bigtiff <- (magic == 43L)
  if (magic != 42L && magic != 43L)
    return(list(descriptions = character(0L), samples_per_pixel = integer(0L)))

  if (bigtiff) {
    readBin(f, "raw", n = 4L)               # bytesize (2) + padding (2)
    ifd_offset <- .read_uint64(f, endian)
  } else {
    ifd_offset <- .ru32(f, endian)
  }

  descriptions      <- character(0L)
  samples_per_pixel <- integer(0L)

  for (pg in seq_len(max_pages)) {
    if (is.null(ifd_offset) || ifd_offset == 0) break

    res <- tryCatch(
      .read_one_ifd_tags(f, ifd_offset, bigtiff, endian),
      error = function(e) NULL
    )
    if (is.null(res)) break

    descriptions      <- c(descriptions,      res$desc %||% "")
    samples_per_pixel <- c(samples_per_pixel, res$samples_per_pixel %||% 1L)
    ifd_offset        <- res$next_offset
  }

  list(descriptions = descriptions, samples_per_pixel = samples_per_pixel)
}

#' Read tags 270 + 277 from one IFD and return the next-IFD offset
#' @keywords internal
.read_one_ifd_tags <- function(f, offset, bigtiff, endian) {
  seek(f, offset)

  n_entries <- if (bigtiff) .read_uint64(f, endian) else
               readBin(f, "integer", n = 1L, size = 2L, signed = FALSE,
                       endian = endian)
  if (is.null(n_entries) || n_entries == 0L) return(NULL)

  entry_size  <- if (bigtiff) 20L else 12L
  vfield_size <- if (bigtiff)  8L else  4L
  entry_start <- seek(f)

  desc      <- NA_character_
  spp       <- 1L
  found_desc <- FALSE; found_spp <- FALSE

  for (e in seq_len(min(n_entries, 500L))) {
    if (found_desc && found_spp) break
    seek(f, entry_start + (e - 1L) * entry_size)
    tag <- readBin(f, "integer", n = 1L, size = 2L, signed = FALSE, endian = endian)
    if (is.null(tag)) break

    if (!tag %in% c(.TAG_IMAGE_DESC, .TAG_SAMPLES_PIX)) next

    type  <- readBin(f, "integer", n = 1L, size = 2L, signed = FALSE, endian = endian)
    count <- if (bigtiff) .read_uint64(f, endian) else {
               v <- readBin(f, "integer", n = 1L, size = 4L, signed = TRUE, endian = endian)
               if (length(v) == 0L) next
               if (v < 0L) as.numeric(v) + 2^32 else as.numeric(v)
             }
    if (is.null(count)) next
    vfield <- readBin(f, "raw", n = vfield_size)

    if (tag == .TAG_IMAGE_DESC && type == 2L) {         # ASCII string
      n_bytes <- count
      if (n_bytes <= vfield_size) {
        raw_str <- vfield[seq_len(n_bytes)]
      } else {
        off <- if (bigtiff) .raw_to_uint64(vfield, endian)
               else         .raw_to_uint32(vfield, endian)
        cur <- seek(f)
        seek(f, off)
        raw_str <- readBin(f, "raw", n = min(n_bytes, 10L * 1024L * 1024L))
        seek(f, cur)
      }
      nuls <- which(raw_str == as.raw(0L))
      if (length(nuls)) raw_str <- raw_str[seq_len(nuls[1L] - 1L)]
      desc <- rawToChar(raw_str)
      found_desc <- TRUE
    } else if (tag == .TAG_SAMPLES_PIX && type == 3L) { # SHORT uint16
      spp <- .raw_to_uint16(vfield[1L:2L], endian)
      found_spp <- TRUE
    }
  }

  # Next IFD offset at end of directory
  seek(f, entry_start + n_entries * entry_size)
  next_offset <- if (bigtiff) .read_uint64(f, endian) else .ru32(f, endian)

  list(desc = desc, samples_per_pixel = spp, next_offset = next_offset)
}

# ============================================================
# Layout detection - channel count and pyramid levels
# ============================================================

#' Determine number of channels and pyramid levels from IFD descriptions
#'
#' Detects the repeating-page pattern:
#' \itemize{
#'   \item All identical -> Polaris; channel count from ScanBands-i.
#'   \item Repeating with period N -> Fusion paged; N channels.
#' }
#' @keywords internal
.detect_layout <- function(descriptions, samples_per_pixel) {
  n <- length(descriptions)
  if (n == 0L) return(list(n_channels = 0L, n_levels = 1L, is_rgb = FALSE))

  is_rgb <- !is.na(samples_per_pixel[1L]) && samples_per_pixel[1L] >= 3L

  # Brightfield: one RGB page per pyramid level; n_ch = samples per pixel.
  # Handle this before XML inspection - RGB pages use a different layout.
  if (is_rgb) {
    spp_val <- as.integer(samples_per_pixel[1L])
    return(list(n_channels = spp_val, n_levels = max(1L, n), is_rgb = TRUE))
  }

  # Strip NA / empty descriptions (non-data pages: label, thumbnail, etc.)
  # NA occurs when a page has no ImageDescription tag (e.g. regular TIFF).
  valid       <- !is.na(descriptions) & nzchar(descriptions)
  clean_descs <- descriptions[valid]

  # No XML at all -> regular non-QPI TIFF; treat each page as one channel
  if (length(clean_descs) == 0L)
    return(list(n_channels = n, n_levels = 1L, is_rgb = is_rgb))

  first <- clean_descs[[1L]]
  m     <- length(clean_descs)

  # All identical -> Polaris (shared XML)
  if (all(clean_descs == first)) {
    n_ch <- tryCatch({
      doc   <- xml2::read_xml(first)
      bands <- xml2::xml_find_all(doc, ".//ScanBands-i")
      if (length(bands) > 0L) length(bands) else 0L
    }, error = function(e) 0L)
    if (n_ch == 0L) {
      # Non-QPI TIFF (e.g. synthetic test files): treat each page as a channel
      return(list(n_channels = n, n_levels = 1L, is_rgb = is_rgb))
    }
    n_lev <- max(1L, m %/% n_ch)
    return(list(n_channels = n_ch, n_levels = n_lev, is_rgb = is_rgb))
  }

  # Find where first description repeats -> period = channel count
  rep_idx <- which(clean_descs == first)
  if (length(rep_idx) < 2L) {
    # Fusion paged: every page has unique per-page XML (no raw-string repetition).
    # Fall back to comparing extracted biomarker names - stop at first repeat.
    n_ch <- .detect_period_from_channel_names(clean_descs)
    if (n_ch > 0L) {
      n_lev <- max(1L, m %/% n_ch)
      return(list(n_channels = n_ch, n_levels = n_lev, is_rgb = is_rgb))
    }
    return(list(n_channels = m, n_levels = 1L, is_rgb = is_rgb))
  }
  period <- rep_idx[2L] - rep_idx[1L]

  # Count full data pages (those whose description matches the cycling pattern)
  n_data <- period
  for (start in seq(period + 1L, m - period + 1L, by = period)) {
    end <- start + period - 1L
    if (end > m) break
    if (identical(clean_descs[start:end], clean_descs[1L:period])) {
      n_data <- end
    } else {
      break
    }
  }

  n_lev <- max(1L, n_data %/% period)
  list(n_channels = period, n_levels = n_lev, is_rgb = is_rgb)
}

# Parse biomarker name from one QPI XML string
.get_qpi_biomarker <- function(xml_str) {
  if (!nzchar(xml_str)) return(NA_character_)
  doc <- tryCatch(xml2::read_xml(xml_str), error = function(e) NULL)
  if (is.null(doc)) return(NA_character_)
  for (tag in .BIOMARKER_TAGS) {
    el <- xml2::xml_find_first(doc, tag)
    if (!inherits(el, "xml_missing")) {
      val <- trimws(xml2::xml_text(el))
      if (nzchar(val) && !val %in% c("None", "--")) return(val)
    }
  }
  NA_character_
}

# Find the repeating period (= n_channels) in Fusion paged XMLs by comparing
# extracted biomarker names rather than raw XML strings.
.detect_period_from_channel_names <- function(descs) {
  if (length(descs) == 0L) return(0L)
  first_name <- .get_qpi_biomarker(descs[[1L]])
  if (is.na(first_name)) return(0L)
  for (i in seq(2L, length(descs))) {
    nm <- .get_qpi_biomarker(descs[[i]])
    if (!is.na(nm) && nm == first_name) return(i - 1L)
  }
  0L
}

# Dimension-based page index computation: groups all IFD pages by (H, W) and
# returns the 1-based page indices that correspond to the requested pyramid
# level and channel subset.  Pages that are thumbnails, labels, or overviews
# (i.e. whose (H, W) group does not contain exactly n_ch pages) are skipped.
.compute_data_page_indices <- function(all_layouts, n_ch, level, ch_idx) {
  hw <- vapply(all_layouts, function(l) {
    h <- as.integer(l$image_height %||% 0L)
    w <- as.integer(l$image_width  %||% 0L)
    paste0(h, "x", w)
  }, character(1L))

  unique_dims <- unique(hw)
  groups      <- lapply(unique_dims, function(k) which(hw == k))
  names(groups) <- unique_dims

  # Keep only groups that could represent a full pyramid level
  data_groups <- groups[vapply(groups, length, integer(1L)) == n_ch]

  if (length(data_groups) == 0L) {
    # No group has exactly n_ch pages; fall back to arithmetic indexing
    return((level - 1L) * n_ch + ch_idx)
  }

  # Sort by decreasing pixel area: full resolution comes first
  areas <- vapply(names(data_groups), function(k) {
    parts <- strsplit(k, "x", fixed = TRUE)[[1L]]
    as.numeric(parts[1L]) * as.numeric(parts[2L])
  }, numeric(1L))
  data_groups <- data_groups[order(areas, decreasing = TRUE)]

  if (level > length(data_groups))
    stop("Requested level ", level, " exceeds available pyramid levels (",
         length(data_groups), ").")

  level_pages <- data_groups[[level]]
  level_pages[ch_idx]
}

# ============================================================
# xml2 helper shortcuts
# ============================================================

# Safe text extraction (direct child)
.xt <- function(node, xpath) {
  el <- xml2::xml_find_first(node, xpath)
  if (inherits(el, "xml_missing")) return(NULL)
  txt <- trimws(xml2::xml_text(el))
  if (!nzchar(txt)) NULL else txt
}
# Float from direct child
.xf <- function(node, xpath) {
  txt <- .xt(node, xpath)
  if (is.null(txt)) return(NULL)
  v <- suppressWarnings(as.numeric(txt))
  if (is.na(v)) NULL else v
}
# Integer from direct child
.xi <- function(node, xpath) {
  v <- .xf(node, xpath)
  if (is.null(v)) NULL else as.integer(v)
}
# Logical from direct child
.xb <- function(node, xpath) {
  txt <- .xt(node, xpath)
  if (is.null(txt)) return(NULL)
  tolower(txt) %in% c("true", "1", "yes")
}
# Float from a node object directly (not xpath)
.xf_node <- function(node) {
  if (is.null(node) || inherits(node, "xml_missing")) return(NULL)
  v <- suppressWarnings(as.numeric(trimws(xml2::xml_text(node))))
  if (is.na(v)) NULL else v
}
# Parse <Color>r,g,b</Color>
.parse_color_node <- function(node) {
  if (is.null(node) || inherits(node, "xml_missing")) return(NULL)
  txt    <- trimws(xml2::xml_text(node))
  parts  <- strsplit(txt, ",", fixed = TRUE)[[1L]]
  if (length(parts) != 3L) return(NULL)
  vals   <- suppressWarnings(as.integer(parts))
  if (anyNA(vals)) NULL else vals
}

# ============================================================
# Low-level binary helpers
# ============================================================

.raw_to_uint16 <- function(raw2, endian) {
  stopifnot(length(raw2) >= 2L)
  b <- as.integer(raw2[1L:2L])
  if (endian == "little") b[1L] + b[2L] * 256L
  else                     b[2L] + b[1L] * 256L
}

.raw_to_uint32 <- function(raw4, endian) {
  stopifnot(length(raw4) >= 4L)
  b <- as.integer(raw4[1L:4L])
  if (endian == "little")
    b[1L] + b[2L] * 256L + b[3L] * 65536L + b[4L] * 16777216L
  else
    b[4L] + b[3L] * 256L + b[2L] * 65536L + b[1L] * 16777216L
}

.raw_to_uint64 <- function(raw8, endian) {
  stopifnot(length(raw8) >= 8L)
  b <- as.numeric(raw8[1L:8L])
  if (endian == "little") sum(b * 256^(0:7))
  else                    sum(b * 256^(7:0))
}

.read_uint64 <- function(f, endian) {
  raw8 <- readBin(f, "raw", n = 8L)
  if (length(raw8) < 8L) return(NULL)
  .raw_to_uint64(raw8, endian)
}

.ru32 <- function(f, endian) {
  v <- readBin(f, "integer", n = 1L, size = 4L, signed = TRUE, endian = endian)
  if (length(v) == 0L) return(NULL)
  # Signed int32 wraps for values > 2^31 - 1; convert to unsigned double
  if (v < 0L) as.numeric(v) + 2^32 else as.numeric(v)
}
