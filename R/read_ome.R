## ============================================================
## read_ome.R
##
## Standard OME-TIFF reader.  Handles files whose first-page ImageDescription
## (tag 270) is an <OME> document (as written by tifffile / Bio-Formats /
## QuPath), where channel names live in the <Channel Name="..."> *attribute*
## rather than the child elements the PerkinElmer QPI parser expects.
##
## Pixel I/O reuses the QPTIFF page readers (.read_qptiff_eager /
## QPTIFFArraySeed): OME-TIFF fluorescence channels are one grayscale TIFF page
## each, identical to the Fusion-paged QPTIFF layout.  The OME-specific work is
## (a) parsing the OME-XML into the shared QPTIFFMetadata model and (b) mapping
## each channel to its absolute IFD page via the <TiffData> table.
## ============================================================

# ---- Detection -----------------------------------------------------------

#' Is a TIFF ImageDescription an OME-XML document?
#' @keywords internal
.is_ome_xml <- function(desc) {
  if (is.null(desc) || length(desc) != 1L || is.na(desc) || !nzchar(desc))
    return(FALSE)
  # Cheap pre-filter before a full parse.
  if (!grepl("<OME", desc, fixed = TRUE) &&
      !grepl("openmicroscopy", desc, fixed = TRUE))
    return(FALSE)
  if (!requireNamespace("xml2", quietly = TRUE)) return(FALSE)
  doc <- tryCatch(xml2::read_xml(desc), error = function(e) NULL)
  if (is.null(doc)) return(FALSE)
  identical(xml2::xml_name(doc), "OME")
}

# ---- OME-XML attribute helpers (namespace stripped) ----------------------

.oa_chr <- function(node, name) {
  if (is.null(node) || inherits(node, "xml_missing")) return(NULL)
  v <- xml2::xml_attr(node, name)
  if (is.na(v) || !nzchar(v)) NULL else v
}
.oa_num <- function(node, name) {
  v <- .oa_chr(node, name)
  if (is.null(v)) return(NULL)
  n <- suppressWarnings(as.numeric(v))
  if (is.na(n)) NULL else n
}
.oa_int <- function(node, name) {
  n <- .oa_num(node, name)
  if (is.null(n)) NULL else as.integer(n)
}

# Parse an OME <Color> attribute (a signed 32-bit RGBA int) into c(r, g, b).
.ome_color_to_rgb <- function(node) {
  v <- .oa_chr(node, "Color")
  if (is.null(v)) return(NULL)
  n <- suppressWarnings(as.numeric(v))
  if (is.na(n)) return(NULL)
  # OME stores RGBA packed big-endian in a signed int32.
  u <- if (n < 0) n + 2^32 else n
  r <- floor(u / 2^24) %% 256
  g <- floor(u / 2^16) %% 256
  b <- floor(u / 2^8)  %% 256
  as.integer(c(r, g, b))
}

# ============================================================
# OME-XML parser
# ============================================================

#' Parse an OME-XML ImageDescription into an intermediate list
#'
#' Namespaces are stripped so plain XPath works.  Reads the first \code{<Image>}
#' (the primary scene): Pixels geometry + physical size, per-channel
#' Name/Fluor/Color/wavelengths, per-Plane exposure/position, the TiffData
#' page map, and the \code{qpi://vectra} MapAnnotation (restored on round-trip).
#'
#' @param xml_str Raw OME-XML string from TIFF tag 270.
#' @return A named list; see \code{.ome_to_metadata} for how it is consumed.
#' @keywords internal
.parse_ome_xml <- function(xml_str) {
  doc <- xml2::read_xml(xml_str)
  xml2::xml_ns_strip(doc)

  image  <- xml2::xml_find_first(doc, ".//Image")
  pixels <- xml2::xml_find_first(doc, ".//Pixels")

  sizes <- list(
    size_x = .oa_int(pixels, "SizeX"), size_y = .oa_int(pixels, "SizeY"),
    size_z = .oa_int(pixels, "SizeZ") %||% 1L,
    size_c = .oa_int(pixels, "SizeC") %||% 1L,
    size_t = .oa_int(pixels, "SizeT") %||% 1L
  )

  # ---- Channels ----------------------------------------------------------
  ch_nodes <- xml2::xml_find_all(pixels, "./Channel")
  channels <- lapply(seq_along(ch_nodes), function(i) {
    n  <- ch_nodes[[i]]
    id <- .oa_chr(n, "ID")                         # "Channel:0:3"
    idx <- suppressWarnings(as.integer(sub(".*:", "", id %||% "")))
    if (is.na(idx)) idx <- i - 1L
    list(
      index                    = idx,
      name                     = .oa_chr(n, "Name"),
      fluorophore              = .oa_chr(n, "Fluor"),
      color_rgb                = .ome_color_to_rgb(n),
      emission_wavelength_nm   = .oa_num(n, "EmissionWavelength"),
      excitation_wavelength_nm = .oa_num(n, "ExcitationWavelength"),
      samples_per_pixel        = .oa_int(n, "SamplesPerPixel") %||% 1L
    )
  })

  # ---- Planes (exposure / stage position keyed by TheC) ------------------
  pl_nodes <- xml2::xml_find_all(pixels, "./Plane")
  planes <- lapply(pl_nodes, function(n) list(
    the_c            = .oa_int(n, "TheC"),
    exposure_time_us = .oa_num(n, "ExposureTime"),
    position_x       = .oa_num(n, "PositionX"),
    position_y       = .oa_num(n, "PositionY")
  ))

  # ---- TiffData page map -------------------------------------------------
  td_nodes <- xml2::xml_find_all(pixels, "./TiffData")
  tiffdata <- lapply(td_nodes, function(n) list(
    first_c     = .oa_int(n, "FirstC") %||% 0L,
    first_z     = .oa_int(n, "FirstZ") %||% 0L,
    first_t     = .oa_int(n, "FirstT") %||% 0L,
    ifd         = .oa_int(n, "IFD")    %||% 0L,
    plane_count = .oa_int(n, "PlaneCount")
  ))

  # ---- qpi://vectra MapAnnotation (round-trip of PerkinElmer fields) -----
  qpi <- list()
  ann_nodes <- xml2::xml_find_all(doc, ".//StructuredAnnotations/MapAnnotation")
  for (a in ann_nodes) {
    if (!identical(.oa_chr(a, "Namespace"), "qpi://vectra")) next
    for (m in xml2::xml_find_all(a, ".//M")) {
      k <- .oa_chr(m, "K"); v <- trimws(xml2::xml_text(m))
      if (!is.null(k) && nzchar(v)) qpi[[k]] <- v
    }
  }
  split       <- .split_qpi_map(qpi)   # image-level vs per-channel (ch<N>_*)
  qpi_image   <- split$image
  qpi_channel <- split$channel

  list(
    sizes            = sizes,
    dimension_order  = .oa_chr(pixels, "DimensionOrder") %||% "XYCZT",
    type             = .oa_chr(pixels, "Type"),
    physical_size_x  = .oa_num(pixels, "PhysicalSizeX"),
    physical_size_y  = .oa_num(pixels, "PhysicalSizeY"),
    physical_unit    = .oa_chr(pixels, "PhysicalSizeXUnit") %||% "\u00b5m",
    image_name       = .oa_chr(image, "Name"),
    acquisition_date = xml2::xml_text(xml2::xml_find_first(doc, ".//Image/AcquisitionDate")),
    channels         = channels,
    planes           = planes,
    tiffdata         = tiffdata,
    qpi_image        = qpi_image,
    qpi_channel      = qpi_channel,
    raw_xml          = xml_str
  )
}

# ---- Channel -> absolute IFD page (1-based) ------------------------------

# Plane index of coordinate (c, z, t) under an OME DimensionOrder string.
.ome_plane_index <- function(c, z, t, sizes, dim_order) {
  dims       <- strsplit(dim_order, "")[[1L]]
  plane_dims <- dims[dims %in% c("C", "Z", "T")]   # last three, fastest first
  coord <- c(C = c, Z = z, T = t)
  size  <- c(C = sizes$size_c, Z = sizes$size_z, T = sizes$size_t)
  idx <- 0L; stride <- 1L
  for (d in plane_dims) {
    idx    <- idx + coord[[d]] * stride
    stride <- stride * size[[d]]
  }
  idx
}

#' Map each channel (Z=T=0) to its 1-based absolute IFD page
#' @keywords internal
.ome_page_indices <- function(ome, n_ch) {
  sizes <- ome$sizes
  dord  <- ome$dimension_order

  # plane index -> IFD (0-based)
  plane_to_ifd <- integer(0L)
  total <- sizes$size_c * sizes$size_z * sizes$size_t
  if (length(ome$tiffdata) > 0L) {
    for (td in ome$tiffdata) {
      p0 <- .ome_plane_index(td$first_c, td$first_z, td$first_t, sizes, dord)
      # OME default for an omitted PlaneCount is "the remaining planes"; when
      # several TiffData entries omit it, each addresses a single plane (the
      # one-TiffData-per-plane convention Bio-Formats / tifffile emit).
      n  <- td$plane_count %||%
              (if (length(ome$tiffdata) == 1L) total - p0 else 1L)
      for (k in seq_len(n) - 1L)
        plane_to_ifd[as.character(p0 + k)] <- td$ifd + k
    }
  }

  vapply(seq_len(n_ch) - 1L, function(ci) {
    p   <- .ome_plane_index(ci, 0L, 0L, sizes, dord)
    # Single-bracket lookup returns NA (not an error) for planes absent from a
    # sparse TiffData table, so the identity fallback below is actually reached.
    ifd <- if (length(plane_to_ifd) > 0L) unname(plane_to_ifd[as.character(p)]) else p
    if (length(ifd) != 1L || is.na(ifd)) ifd <- p   # default identity mapping
    as.integer(ifd) + 1L
  }, integer(1L))
}

# ============================================================
# OME-XML -> QPTIFFMetadata
# ============================================================

# Coerce a named list of string qpi fields to typed values for a channel.
.coerce_qpi_channel <- function(qc) {
  num <- c("responsivity")
  int <- c("bit_depth", "offset_counts", "signal_units",
           "roi_x", "roi_y", "roi_width", "roi_height")
  lgl <- c("is_unmixed_component", "autofluorescence_subtracted")
  out <- list()
  for (k in names(qc)) {
    v <- qc[[k]]
    out[[k]] <- if (k %in% num) suppressWarnings(as.numeric(v))
                else if (k %in% int) suppressWarnings(as.integer(as.numeric(v)))
                else if (k %in% lgl) tolower(v) %in% c("true", "1", "yes")
                else v
  }
  out
}

#' Build a QPTIFFMetadata object from a parsed OME-XML list
#' @keywords internal
.ome_to_metadata <- function(ome, n_levels = 1L, format = .QPTIFF_OME_TIFF) {
  qi <- ome$qpi_image

  slide <- .new_slide_info(
    slide_id             = qi$slide_id,
    barcode              = qi$barcode,
    study_name           = qi$study_name,
    computer_name        = qi$computer_name,
    datetime             = qi$datetime %||%
                             (if (nzchar(ome$acquisition_date %||% "")) ome$acquisition_date else NULL),
    acquisition_software = qi$acquisition_software,
    description_version  = qi$description_version,
    identifier           = qi$identifier
  )

  camera <- .new_camera_info(
    camera_name = qi$camera_name,
    gain        = if (!is.null(qi$camera_gain)) suppressWarnings(as.numeric(qi$camera_gain)) else NULL,
    bit_depth   = if (!is.null(qi$camera_bit_depth)) suppressWarnings(as.integer(qi$camera_bit_depth)) else NULL
  )
  scan_res <- .new_scan_resolution(base_pixel_size_um = ome$physical_size_x)
  image_info <- .new_image_info(
    image_type        = qi$image_type,
    objective         = qi$objective,
    bf_lamp_type      = qi$bf_lamp_type,
    scan_profile_name = qi$scan_profile_name,
    scan_mode         = qi$scan_mode,
    camera            = camera,
    scan_resolution   = scan_res
  )
  # Preserve the original acquisition format + any unmapped vendor keys.
  keep <- setdiff(names(qi),
                  c("slide_id", "barcode", "study_name", "computer_name",
                    "datetime", "acquisition_software", "description_version",
                    "identifier", "camera_name", "camera_gain",
                    "camera_bit_depth", "image_type", "objective",
                    "bf_lamp_type", "scan_profile_name", "scan_mode"))
  if (length(keep) > 0L) image_info$extra <- qi[keep]

  planes_by_c <- list()
  for (pl in ome$planes)
    if (!is.null(pl$the_c)) planes_by_c[[as.character(pl$the_c)]] <- pl

  channels <- lapply(ome$channels, function(ch) {
    out <- list(
      index                    = ch$index,
      name                     = ch$name,
      fluorophore              = ch$fluorophore,
      color_rgb                = ch$color_rgb,
      emission_wavelength_nm   = ch$emission_wavelength_nm,
      excitation_wavelength_nm = ch$excitation_wavelength_nm,
      is_brightfield           = FALSE
    )
    pl <- planes_by_c[[as.character(ch$index)]]
    if (!is.null(pl) && !is.null(pl$exposure_time_us))
      out$exposure_time_us <- pl$exposure_time_us
    qc <- ome$qpi_channel[[as.character(ch$index)]]
    if (!is.null(qc)) out <- utils::modifyList(out, .coerce_qpi_channel(qc))
    # Split a compound biomarker name (CD3e-AF647) into protein + dye.
    if (!is.null(out$name)) {
      sp <- .split_protein_dye(out$name, out$fluorophore)
      out$name <- sp$protein
      out$dye_from_name <- sp$dye
    }
    out
  })

  .single_scene_metadata(
    slide      = slide,
    image_info = image_info,
    channels   = channels,
    scales     = .build_scales(n_levels, ome$physical_size_x),
    format     = format,
    raw_xml    = ome$raw_xml
  )
}

# ============================================================
# Top-level OME-TIFF reader
# ============================================================

#' Read a standard OME-TIFF into a QPTIFFImage
#'
#' @param path      File path.
#' @param ifd_info  Result of \code{.read_all_ifd_info(path)} (already read by
#'   \code{read_qptiff} during dispatch).
#' @inheritParams read_qptiff
#' @keywords internal
.read_ome_tiff <- function(path, ifd_info, channels = NULL, level = 1L,
                           as_integer = TRUE, lazy = FALSE) {
  if (!requireNamespace("xml2", quietly = TRUE))
    stop("Package 'xml2' is required to read OME-TIFF files.")
  if (level != 1L)
    stop("OME-TIFF sub-resolution pyramids (SubIFDs) are not supported; ",
         "use level = 1.")

  message("Reading OME-TIFF metadata ...")
  ome  <- .parse_ome_xml(ifd_info$descriptions[[1L]])
  n_pg <- length(ifd_info$descriptions)
  n_ch <- ome$sizes$size_c %||% n_pg
  if (is.null(n_ch) || n_ch < 1L) n_ch <- n_pg
  n_ch <- min(as.integer(n_ch), n_pg)

  # Channel names: <Channel Name> if present, else Channel_i.
  all_channel_names <- vapply(ome$channels,
                              function(c) c$name %||% NA_character_,
                              character(1L))
  if (length(all_channel_names) < n_ch)
    all_channel_names <- c(all_channel_names,
                           rep(NA_character_, n_ch - length(all_channel_names)))
  all_channel_names <- all_channel_names[seq_len(n_ch)]
  miss <- is.na(all_channel_names) | !nzchar(all_channel_names)
  all_channel_names[miss] <- paste0("Channel_", which(miss) - 1L)

  page_of_channel <- .ome_page_indices(ome, n_ch)

  # --- Select channels ---------------------------------------------------
  sel          <- .select_channels(channels, all_channel_names)
  ch_idx       <- sel$ch_idx
  ch_names     <- sel$ch_names
  page_indices <- page_of_channel[ch_idx]

  message("Reading IFD page layouts ...")
  all_layouts_raw <- .read_all_ifd_page_layouts(path)

  meta <- .ome_to_metadata(ome, n_levels = 1L, format = .QPTIFF_OME_TIFF)
  # Keep the stored channel metadata aligned 1:1 with the loaded channels.
  meta <- .subset_scene_channels(meta, all_channel_names, ch_idx)

  message("Loading ", length(ch_idx), " channel(s) ...")
  if (lazy)
    .read_qptiff_lazy(path, ch_idx, ch_names, page_indices,
                      all_layouts_raw, level, meta, as_integer)
  else
    .read_qptiff_eager(path, ch_idx, ch_names, page_indices,
                       all_layouts_raw, meta, as_integer, is_rgb = FALSE)
}
