## Tests for the rewritten QPTIFF reader (ported from Python bioio-tifffile)

test_that("read_qptiff throws error for missing file", {
  expect_error(read_qptiff("nonexistent.qptiff"), "File not found")
})

# ---------- Binary helper tests ------------------------------------------

test_that(".raw_to_uint16 converts correctly (little-endian)", {
  # 300 = 0x012C  in little-endian: 2C 01
  raw2 <- as.raw(c(0x2C, 0x01))
  expect_equal(bgnormR:::.raw_to_uint16(raw2, "little"), 300L)
})

test_that(".raw_to_uint16 converts correctly (big-endian)", {
  # 300 = 0x012C  in big-endian: 01 2C
  raw2 <- as.raw(c(0x01, 0x2C))
  expect_equal(bgnormR:::.raw_to_uint16(raw2, "big"), 300L)
})

test_that(".raw_to_uint32 converts correctly (little-endian)", {
  # 256 = 0x00000100 in little-endian: 00 01 00 00
  raw4 <- as.raw(c(0x00, 0x01, 0x00, 0x00))
  expect_equal(bgnormR:::.raw_to_uint32(raw4, "little"), 256L)
})

test_that(".raw_to_uint32 converts correctly (big-endian)", {
  # 256 = 0x00000100 in big-endian: 00 00 01 00
  raw4 <- as.raw(c(0x00, 0x00, 0x01, 0x00))
  expect_equal(bgnormR:::.raw_to_uint32(raw4, "big"), 256L)
})

test_that(".raw_to_uint64 handles values above 2^32", {
  # 2^33 = 8589934592  little-endian: 00 00 00 00 02 00 00 00
  raw8 <- as.raw(c(0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00))
  expect_equal(bgnormR:::.raw_to_uint64(raw8, "little"), 2^33)
})

# ---------- Format detection tests ---------------------------------------

test_that(".detect_qptiff_format returns brightfield for Brightfield scan_mode", {
  skip_if_not_installed("xml2")
  doc <- xml2::read_xml("<PerkinElmer-QPI-ImageDescription/>")
  expect_equal(bgnormR:::.detect_qptiff_format(doc, "Brightfield", FALSE),
               "brightfield")
})

test_that(".detect_qptiff_format returns brightfield for RGB pages", {
  skip_if_not_installed("xml2")
  doc <- xml2::read_xml("<PerkinElmer-QPI-ImageDescription/>")
  expect_equal(bgnormR:::.detect_qptiff_format(doc, "fluorescence", TRUE),
               "brightfield")
})

test_that(".detect_qptiff_format returns polaris_scanband for ScanBands-i XML", {
  skip_if_not_installed("xml2")
  xml <- '<PerkinElmer-QPI-ImageDescription>
    <ScanProfile>
      <root>
        <ScanBands-i><Biomarker>DAPI</Biomarker></ScanBands-i>
        <ScanBands-i><Biomarker>CD3</Biomarker></ScanBands-i>
      </root>
    </ScanProfile>
  </PerkinElmer-QPI-ImageDescription>'
  doc <- xml2::read_xml(xml)
  expect_equal(bgnormR:::.detect_qptiff_format(doc, "fluorescence", FALSE),
               "polaris_scanband")
})

test_that(".detect_qptiff_format returns fusion_paged for JSON ScanProfile", {
  skip_if_not_installed("xml2")
  xml <- '<PerkinElmer-QPI-ImageDescription>
    <ScanProfile>{"isTma": false, "binning": 1}</ScanProfile>
  </PerkinElmer-QPI-ImageDescription>'
  doc <- xml2::read_xml(xml)
  expect_equal(bgnormR:::.detect_qptiff_format(doc, "fluorescence", FALSE),
               "fusion_paged")
})

# ---------- Channel parsing tests ----------------------------------------

test_that(".parse_qpi_xml extracts channels from Polaris ScanBands-i XML", {
  skip_if_not_installed("xml2")
  polaris_xml <- '<?xml version="1.0"?>
  <PerkinElmer-QPI-ImageDescription>
    <ScanProfile>
      <root>
        <ScanBands-i>
          <Biomarker>DAPI</Biomarker>
          <Fluorophore>DAPI</Fluorophore>
          <ExposureTime>30000</ExposureTime>
        </ScanBands-i>
        <ScanBands-i>
          <Biomarker>CD3</Biomarker>
          <Fluorophore>OPAL520</Fluorophore>
          <ExposureTime>200000</ExposureTime>
        </ScanBands-i>
        <ScanBands-i>
          <Biomarker>CD8</Biomarker>
          <Fluorophore>OPAL650</Fluorophore>
          <ExposureTime>150000</ExposureTime>
        </ScanBands-i>
      </root>
    </ScanProfile>
  </PerkinElmer-QPI-ImageDescription>'
  meta <- bgnormR:::.parse_qpi_xml(polaris_xml, per_page_xmls = NULL,
                                    n_channels = 3L, is_rgb = FALSE)
  expect_equal(meta$format, "polaris_scanband")
  expect_length(meta$channels, 3L)
  expect_equal(meta$channels[[1]]$name, "DAPI")
  expect_equal(meta$channels[[2]]$name, "CD3")
  expect_equal(meta$channels[[3]]$name, "CD8")
  expect_equal(meta$channels[[2]]$fluorophore, "OPAL520")
  expect_equal(meta$channels[[1]]$exposure_time_us, 30000)
})

test_that(".parse_qpi_xml extracts channels from Fusion per-page XMLs", {
  skip_if_not_installed("xml2")
  root_xml <- '<?xml version="1.0"?>
  <PerkinElmer-QPI-ImageDescription>
    <ScanProfile>{"isTma": false, "binning": 1}</ScanProfile>
  </PerkinElmer-QPI-ImageDescription>'

  make_page_xml <- function(biomarker, fluor, et) {
    sprintf('<PerkinElmer-QPI-ImageDescription>
      <Biomarker>%s</Biomarker>
      <Fluorophore>%s</Fluorophore>
      <ExposureTime>%d</ExposureTime>
      <ImageType>FullResolution</ImageType>
    </PerkinElmer-QPI-ImageDescription>', biomarker, fluor, et)
  }
  pp_xmls <- c(
    make_page_xml("DAPI", "DAPI",   30000L),
    make_page_xml("CD3",  "AF488", 200000L),
    make_page_xml("CD8",  "AF555", 150000L)
  )

  meta <- bgnormR:::.parse_qpi_xml(root_xml, per_page_xmls = pp_xmls,
                                    n_channels = 3L, is_rgb = FALSE)
  expect_equal(meta$format, "fusion_paged")
  expect_length(meta$channels, 3L)
  expect_equal(meta$channels[[1]]$name, "DAPI")
  expect_equal(meta$channels[[3]]$name, "CD8")
  expect_equal(meta$channels[[2]]$fluorophore, "AF488")
  expect_equal(meta$channels[[1]]$exposure_time_us, 30000)
})

test_that(".parse_qpi_xml handles brightfield XML", {
  skip_if_not_installed("xml2")
  bf_xml <- '<?xml version="1.0"?>
  <PerkinElmer-QPI-ImageDescription>
    <ImageType>FullResolution</ImageType>
    <ScanProfile>
      <root>
        <Mode>Brightfield</Mode>
      </root>
    </ScanProfile>
  </PerkinElmer-QPI-ImageDescription>'
  meta <- bgnormR:::.parse_qpi_xml(bf_xml, n_channels = 3L, is_rgb = TRUE)
  expect_equal(meta$format, "brightfield")
  expect_length(meta$channels, 3L)
  expect_equal(meta$channels[[1]]$name, "Red")
  expect_true(meta$channels[[1]]$is_brightfield)
})

test_that(".populate_channel_from_node extracts wavelength from filter bands", {
  skip_if_not_installed("xml2")
  node_xml <- '<Band>
    <Biomarker>CD3</Biomarker>
    <Fluorophore>OPAL520</Fluorophore>
    <EmissionFilter>
      <Bands>
        <Band>
          <Active>true</Active>
          <Cuton>495</Cuton>
          <Cutoff>545</Cutoff>
        </Band>
      </Bands>
    </EmissionFilter>
    <ExcitationFilter>
      <Bands>
        <Band>
          <Active>true</Active>
          <Cuton>465</Cuton>
          <Cutoff>495</Cutoff>
        </Band>
      </Bands>
    </ExcitationFilter>
  </Band>'
  node <- xml2::read_xml(node_xml)
  fields <- list(index = 0L, name = "Channel_0", is_brightfield = FALSE)
  out <- bgnormR:::.populate_channel_from_node(node, fields)
  expect_equal(out$name, "CD3")
  expect_equal(out$fluorophore, "OPAL520")
  expect_equal(out$emission_wavelength_nm, (495 + 545) / 2)
  expect_equal(out$excitation_wavelength_nm, (465 + 495) / 2)
})

test_that(".populate_channel_from_node infers wavelength from OPAL fluorophore name", {
  skip_if_not_installed("xml2")
  node_xml <- '<Band>
    <Biomarker>CD8</Biomarker>
    <Fluorophore>OPAL650</Fluorophore>
  </Band>'
  node <- xml2::read_xml(node_xml)
  fields <- list(index = 1L, name = "Channel_1", is_brightfield = FALSE)
  out <- bgnormR:::.populate_channel_from_node(node, fields)
  expect_equal(out$emission_wavelength_nm, 650)
})

# ---------- Protein / dye splitting tests --------------------------------

test_that(".split_protein_dye handles simple Alexa Fluor suffix", {
  res <- bgnormR:::.split_protein_dye("CD209-AF647")
  expect_equal(res$protein, "CD209")
  expect_equal(res$dye,     "AF647")
})

test_that(".split_protein_dye handles Atto suffix", {
  res <- bgnormR:::.split_protein_dye("CD31-Atto550")
  expect_equal(res$protein, "CD31")
  expect_equal(res$dye,     "Atto550")
})

test_that(".split_protein_dye handles multi-hyphen protein names (HLA-A)", {
  res <- bgnormR:::.split_protein_dye("HLA-A-Atto550")
  expect_equal(res$protein, "HLA-A")
  expect_equal(res$dye,     "Atto550")
})

test_that(".split_protein_dye uses fluorophore field for primary split", {
  # fluorophore matches the suffix exactly → reliable primary path
  res <- bgnormR:::.split_protein_dye("CD3e-AF647", fluorophore = "AF647")
  expect_equal(res$protein, "CD3e")
  expect_equal(res$dye,     "AF647")
})

test_that(".split_protein_dye ignores mismatched fluorophore and falls back to regex", {
  # fluorophore tag says 'CY5' (detection channel) but name suffix is 'AF750'
  res <- bgnormR:::.split_protein_dye("CD107a-AF750", fluorophore = "CY5")
  expect_equal(res$protein, "CD107a")
  expect_equal(res$dye,     "AF750")
})

test_that(".split_protein_dye returns NULL dye for bare protein names", {
  res <- bgnormR:::.split_protein_dye("DAPI")
  expect_equal(res$protein, "DAPI")
  expect_null(res$dye)
})

test_that(".split_protein_dye handles OPAL suffix (Polaris)", {
  res <- bgnormR:::.split_protein_dye("CD8-OPAL650")
  expect_equal(res$protein, "CD8")
  expect_equal(res$dye,     "OPAL650")
})

test_that(".split_protein_dye handles names with spaces", {
  res <- bgnormR:::.split_protein_dye("Histone H3 Phospho Ser28-AF647")
  expect_equal(res$protein, "Histone H3 Phospho Ser28")
  expect_equal(res$dye,     "AF647")
})

test_that(".split_protein_dye handles BLANK channels", {
  res <- bgnormR:::.split_protein_dye("B3-BLANK-AF647")
  expect_equal(res$protein, "B3-BLANK")
  expect_equal(res$dye,     "AF647")
})

test_that("read_qptiff strips dye from channel names in the returned object", {
  skip_if_not_installed("xml2")
  root_xml <- '<?xml version="1.0"?>
  <PerkinElmer-QPI-ImageDescription>
    <ScanProfile>{"isTma": false}</ScanProfile>
  </PerkinElmer-QPI-ImageDescription>'
  make_xml <- function(bm, fl) sprintf(
    '<PerkinElmer-QPI-ImageDescription>
       <Biomarker>%s</Biomarker><Fluorophore>%s</Fluorophore>
     </PerkinElmer-QPI-ImageDescription>', bm, fl)
  pp <- c(make_xml("DAPI", "DAPI"),
          make_xml("CD3e-AF647", "AF647"),
          make_xml("HLA-A-Atto550", "ATTO550"))
  meta <- bgnormR:::.parse_qpi_xml(root_xml, per_page_xmls = pp, n_channels = 3L,
                                    is_rgb = FALSE)
  # Manually apply splitting (as read_qptiff() does)
  meta$channels <- lapply(meta$channels, function(ch) {
    sp <- bgnormR:::.split_protein_dye(ch$name, ch$fluorophore)
    ch$name <- sp$protein; ch$dye_from_name <- sp$dye; ch
  })
  nms <- vapply(meta$channels, `[[`, character(1L), "name")
  dye <- vapply(meta$channels, function(ch) ch$dye_from_name %||% "", character(1L))
  expect_equal(nms, c("DAPI", "CD3e", "HLA-A"))
  expect_equal(dye, c("", "AF647", "Atto550"))
})

test_that("metadata() returns the metadata list from a QPTIFFImage", {
  arr  <- array(0, dim = c(2, 2, 1),
                dimnames = list(NULL, NULL, "DAPI"))
  img  <- bgnormR:::.new_QPTIFFImage(arr, metadata = list(format = "fusion_paged",
                                                            channels = list()))
  meta <- metadata(img)
  expect_equal(meta$format, "fusion_paged")
})

# ---------- Layout detection tests ---------------------------------------

test_that(".detect_layout finds n_channels from ScanBands-i (all-same Polaris)", {
  skip_if_not_installed("xml2")
  polaris_xml <- '<PerkinElmer-QPI-ImageDescription>
    <ScanProfile>
      <root>
        <ScanBands-i><Biomarker>DAPI</Biomarker></ScanBands-i>
        <ScanBands-i><Biomarker>CD3</Biomarker></ScanBands-i>
      </root>
    </ScanProfile>
  </PerkinElmer-QPI-ImageDescription>'
  # Simulate 3 levels × 2 channels = 6 identical descriptions
  descs <- rep(polaris_xml, 6L)
  layout <- bgnormR:::.detect_layout(descs, rep(1L, 6L))
  expect_equal(layout$n_channels, 2L)
  expect_equal(layout$n_levels, 3L)
})

test_that(".detect_layout finds period from Fusion paged format", {
  skip_if_not_installed("xml2")
  # 3 distinct per-page XMLs, repeated over 2 levels → 6 IFDs
  make_xml <- function(bm) sprintf(
    '<PerkinElmer-QPI-ImageDescription><Biomarker>%s</Biomarker></PerkinElmer-QPI-ImageDescription>',
    bm)
  descs <- rep(c(make_xml("DAPI"), make_xml("CD3"), make_xml("CD8")), 2L)
  layout <- bgnormR:::.detect_layout(descs, rep(1L, 6L))
  expect_equal(layout$n_channels, 3L)
  expect_equal(layout$n_levels, 2L)
})

test_that(".detect_layout detects RGB from samples_per_pixel", {
  skip_if_not_installed("xml2")
  descs <- rep("<PerkinElmer-QPI-ImageDescription/>", 3L)
  layout <- bgnormR:::.detect_layout(descs, rep(3L, 3L))
  expect_true(layout$is_rgb)
})

# ---------- QPTIFFImage S3 class tests -----------------------------------

test_that(".new_QPTIFFImage creates a 3-D array with correct class", {
  arr <- array(1:24, dim = c(3, 4, 2),
               dimnames = list(NULL, NULL, c("DAPI", "CD3")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  expect_s3_class(img, "QPTIFFImage")
  expect_true(is.array(img))
  expect_equal(dim(img), c(3L, 4L, 2L))
})

test_that("dim() on QPTIFFImage returns c(H, W, C)", {
  arr <- array(0, dim = c(2, 3, 3),
               dimnames = list(NULL, NULL, c("A", "B", "C")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  expect_equal(dim(img), c(2L, 3L, 3L))
})

test_that("[.QPTIFFImage returns a 1-channel QPTIFFImage for single channel", {
  arr <- array(seq_len(8), dim = c(2, 2, 2),
               dimnames = list(NULL, NULL, c("X", "Y")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  ch  <- img[, , "X"]
  expect_s3_class(ch, "QPTIFFImage")
  expect_equal(dim(ch), c(2L, 2L, 1L))
  expect_equal(dimnames(ch)[[3L]], "X")
  expect_equal(as.vector(ch), as.vector(arr[, , 1L]))
})

test_that("[.QPTIFFImage spatial crop returns QPTIFFImage with correct dims", {
  arr <- array(0, dim = c(10, 10, 3),
               dimnames = list(NULL, NULL, c("A", "B", "C")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  crop <- img[1:5, 1:6, ]
  expect_s3_class(crop, "QPTIFFImage")
  expect_equal(dim(crop), c(5L, 6L, 3L))
})

test_that("[.QPTIFFImage spatial + channel crop returns QPTIFFImage", {
  arr <- array(seq_len(3 * 4 * 2), dim = c(3, 4, 2),
               dimnames = list(NULL, NULL, c("DAPI", "CD3")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  sub <- img[1:2, 1:3, "DAPI"]
  expect_s3_class(sub, "QPTIFFImage")
  expect_equal(dim(sub), c(2L, 3L, 1L))
})

test_that("[,, on QPTIFFImage preserves class for multi-channel subset", {
  arr <- array(0, dim = c(2, 2, 3),
               dimnames = list(NULL, NULL, c("A", "B", "C")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  sub <- img[, , c("A", "C"), drop = FALSE]
  expect_s3_class(sub, "QPTIFFImage")
  expect_equal(dim(sub), c(2L, 2L, 2L))
  expect_equal(dimnames(sub)[[3L]], c("A", "C"))
})

test_that("as.list() returns named list of 2-D matrices", {
  arr <- array(seq_len(12), dim = c(2, 3, 2),
               dimnames = list(NULL, NULL, c("DAPI", "CD3")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  lst <- as.list(img)
  expect_false(inherits(lst, "QPTIFFImage"))
  expect_named(lst, c("DAPI", "CD3"))
  expect_true(is.matrix(lst$DAPI))
  expect_equal(lst$DAPI, arr[, , 1L])
})

test_that("names() and length() return channel names and count", {
  arr <- array(0, dim = c(2, 2, 2),
               dimnames = list(NULL, NULL, c("DAPI", "CD8")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  expect_equal(names(img), c("DAPI", "CD8"))
  expect_equal(length(img), 2L)
})

test_that("names<- replaces channel names and preserves class and metadata", {
  arr <- array(0, dim = c(2, 2, 2),
               dimnames = list(NULL, NULL, c("ch1", "ch2")))
  img <- bgnormR:::.new_QPTIFFImage(arr, metadata = list(foo = "bar"))
  names(img) <- c("DAPI", "CD3")
  expect_s3_class(img, "QPTIFFImage")
  expect_equal(names(img), c("DAPI", "CD3"))
  expect_equal(attr(img, "metadata")$foo, "bar")
})

test_that("names<- errors when length does not match channel count", {
  arr <- array(0, dim = c(2, 2, 3),
               dimnames = list(NULL, NULL, c("a", "b", "c")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  expect_error(names(img) <- c("X", "Y"), "channels")
})

test_that("dimnames<- replaces all dimnames and preserves class and metadata", {
  arr <- array(0, dim = c(2, 3, 2),
               dimnames = list(NULL, NULL, c("A", "B")))
  img <- bgnormR:::.new_QPTIFFImage(arr, metadata = list(x = 1L))
  dimnames(img) <- list(paste0("r", 1:2), paste0("c", 1:3), c("DAPI", "CD8"))
  expect_s3_class(img, "QPTIFFImage")
  expect_equal(names(img), c("DAPI", "CD8"))
  expect_equal(dimnames(img)[[1L]], c("r1", "r2"))
  expect_equal(attr(img, "metadata")$x, 1L)
})

test_that("as.matrix() converts a single-channel QPTIFFImage to a matrix", {
  arr <- array(seq_len(6), dim = c(2, 3, 1),
               dimnames = list(NULL, NULL, "DAPI"))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  m   <- as.matrix(img)
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 3L))
  expect_equal(as.vector(m), seq_len(6))
  expect_equal(attr(m, "channel"), "DAPI")
})

test_that("as.matrix() errors for multi-channel QPTIFFImage", {
  arr <- array(0, dim = c(2, 2, 2),
               dimnames = list(NULL, NULL, c("DAPI", "CD3")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  expect_error(as.matrix(img), "single-channel")
})

test_that("print.QPTIFFImage shows dimensions and channel names", {
  arr <- array(seq_len(60), dim = c(5, 6, 2),
               dimnames = list(NULL, NULL, c("DAPI", "CD3")))
  img <- bgnormR:::.new_QPTIFFImage(arr)
  out <- capture.output(print(img))
  expect_true(any(grepl("QPTIFFImage", out)))
  expect_true(any(grepl("5", out)))
  expect_true(any(grepl("DAPI", out)))
})

# ---------- write_qptiff tests -----------------------------------------------

test_that("write_qptiff writes a readable 16-bit TIFF for a raw QPTIFFImage", {
  skip_if_not_installed("tiff")
  arr <- array(as.double(0:179), dim = c(10, 6, 3),
               dimnames = list(NULL, NULL, c("DAPI", "CD3", "CD8")))
  img <- as.QPTIFFImage(arr)
  tmp <- tempfile(fileext = ".tif")
  on.exit(unlink(tmp), add = TRUE)

  out_path <- write_qptiff(img, tmp)
  expect_equal(out_path, tmp)
  expect_true(file.exists(tmp))

  pages <- tiff::readTIFF(tmp, all = TRUE, as.is = TRUE)
  expect_length(pages, 3L)
  expect_equal(dim(pages[[1L]]), c(10L, 6L))
  expect_equal(max(abs(pages[[1L]] - arr[, , 1L])), 0)
})

test_that("write_qptiff applies 2^x for bgnorm-adjusted QPTIFFImage", {
  skip_if_not_installed("tiff")
  img <- sim_pixel_image(n_markers = 2L, ch_names = c("DAPI", "CD3"))
  res <- bgnorm_pixels(img)
  tmp <- tempfile(fileext = ".tif")
  on.exit(unlink(tmp), add = TRUE)

  write_qptiff(res, tmp)
  pages <- tiff::readTIFF(tmp, all = TRUE, as.is = TRUE)
  expect_length(pages, 2L)

  # Values should equal round(2^x) for background-adjusted log intensities
  log_ch1 <- as.array(res)[, , 1L]
  expected <- pmin(pmax(round(2^log_ch1), 0), 65535)
  expect_equal(max(abs(pages[[1L]] - expected)), 0)
})

test_that("write_qptiff embeds bgnorm metadata as JSON in ImageDescription", {
  img <- sim_pixel_image()
  res <- bgnorm_pixels(img)
  tmp <- tempfile(fileext = ".tif")
  on.exit(unlink(tmp), add = TRUE)

  write_qptiff(res, tmp)
  ifd_info <- bgnormR:::.read_all_ifd_info(tmp)
  desc <- ifd_info$descriptions[[1L]]

  expect_true(nzchar(desc))
  expect_true(grepl("bgnorm", desc))
  expect_true(grepl("jsd",    desc))
  expect_true(grepl("Ch1",    desc))
  expect_true(grepl("2\\^x",  desc))
})

test_that("write_qptiff writes correct channel names in metadata", {
  arr <- array(0, dim = c(4, 4, 2),
               dimnames = list(NULL, NULL, c("PanCK", "CD45")))
  img <- as.QPTIFFImage(arr)
  tmp <- tempfile(fileext = ".tif")
  on.exit(unlink(tmp), add = TRUE)

  write_qptiff(img, tmp)
  ifd_info <- bgnormR:::.read_all_ifd_info(tmp)
  expect_true(grepl("PanCK", ifd_info$descriptions[[1L]]))
  expect_true(grepl("CD45",  ifd_info$descriptions[[2L]]))
})

test_that("write_qptiff errors for non-QPTIFFImage input", {
  expect_error(write_qptiff(matrix(1:4, 2, 2), "out.tif"), "QPTIFFImage")
})

test_that("write_qptiff errors for invalid path", {
  img <- as.QPTIFFImage(array(0, dim = c(2, 2, 1),
                               dimnames = list(NULL, NULL, "Ch1")))
  expect_error(write_qptiff(img, ""), "non-empty")
})

# ---------- QPTIFFMetadata model + accessors -----------------------------

test_that("QPTIFFMetadata accessors read the OME hierarchy", {
  ch <- list(
    list(index = 0L, name = "DAPI", fluorophore = "DAPI",
         exposure_time_us = 3000, color_rgb = c(0L, 0L, 255L)),
    list(index = 1L, name = "CD3", fluorophore = "OPAL520",
         exposure_time_us = 150000, is_brightfield = FALSE)
  )
  meta <- bgnormR:::.single_scene_metadata(
    slide      = bgnormR:::.new_slide_info(slide_id = "S1"),
    image_info = bgnormR:::.new_image_info(
      scan_resolution = bgnormR:::.new_scan_resolution(base_pixel_size_um = 0.5)),
    channels   = ch,
    scales     = bgnormR:::.build_scales(3L, 0.5),
    format     = "fusion_paged")

  expect_s3_class(meta, "QPTIFFMetadata")
  expect_equal(qpi_format(meta), "fusion_paged")
  expect_equal(qpi_channel_names(meta), c("DAPI", "CD3"))
  expect_equal(qpi_pixel_size_um(meta), 0.5)
  expect_equal(qpi_n_levels(meta), 3L)
  ct <- channel_table(meta)
  expect_equal(nrow(ct), 2L)
  expect_equal(ct$name, c("DAPI", "CD3"))
  expect_equal(ct$color_rgb[1L], "0,0,255")
})

# ---------- OME colour encode / decode -----------------------------------

test_that("OME colour round-trips rgb -> int32 -> rgb", {
  skip_if_not_installed("xml2")
  rgb <- c(18L, 52L, 200L)
  enc <- bgnormR:::.rgb_to_ome_color(rgb)
  node <- xml2::read_xml(sprintf('<Channel Color="%s"/>', enc))
  expect_equal(bgnormR:::.ome_color_to_rgb(node), rgb)
})

# ---------- OME-TIFF plane index / detection -----------------------------

test_that(".ome_plane_index follows DimensionOrder", {
  sizes <- list(size_c = 4L, size_z = 1L, size_t = 1L)
  idx <- vapply(0:3, function(c)
    bgnormR:::.ome_plane_index(c, 0L, 0L, sizes, "XYCZT"), integer(1L))
  expect_equal(idx, 0:3)
})

test_that(".is_ome_xml detects an OME root but not a QPI root", {
  skip_if_not_installed("xml2")
  expect_true(bgnormR:::.is_ome_xml(
    '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"/>'))
  expect_false(bgnormR:::.is_ome_xml("<PerkinElmer-QPI-ImageDescription/>"))
  expect_false(bgnormR:::.is_ome_xml(NA_character_))
})

# ---------- OME-TIFF reading (channel Name attributes) -------------------

test_that("read_qptiff reads OME-TIFF channel Name attributes written by write_qptiff", {
  skip_if_not_installed("tiff")
  skip_if_not_installed("xml2")
  arr <- array(as.double(seq_len(4 * 5 * 3)), dim = c(4, 5, 3),
               dimnames = list(NULL, NULL, c("DAPI", "CD3", "CD8")))
  img <- as.QPTIFFImage(arr)
  tmp <- tempfile(fileext = ".ome.tif")
  on.exit(unlink(tmp), add = TRUE)
  write_qptiff(img, tmp)

  img2 <- read_qptiff(tmp)
  expect_equal(qpi_format(metadata(img2)), "ome_tiff")
  expect_equal(names(img2), c("DAPI", "CD3", "CD8"))
  expect_equal(dim(img2), c(4L, 5L, 3L))
  expect_equal(as.numeric(as.array(img2)), as.numeric(arr))
})

test_that("write/read round-trips rich channel + slide metadata via OME + qpi map", {
  skip_if_not_installed("tiff")
  skip_if_not_installed("xml2")
  ch <- list(
    list(index = 0L, name = "DAPI", fluorophore = "DAPI",
         color_rgb = c(0L, 0L, 255L), exposure_time_us = 3000,
         signal_units = 64L, is_brightfield = FALSE),
    list(index = 1L, name = "CD3", fluorophore = "OPAL520",
         color_rgb = c(0L, 255L, 0L), exposure_time_us = 150000,
         emission_wavelength_nm = 520, is_brightfield = FALSE)
  )
  meta <- bgnormR:::.single_scene_metadata(
    slide      = bgnormR:::.new_slide_info(slide_id = "SLIDE1"),
    image_info = bgnormR:::.new_image_info(
      scan_resolution = bgnormR:::.new_scan_resolution(base_pixel_size_um = 0.5)),
    channels   = ch,
    scales     = bgnormR:::.build_scales(1L, 0.5),
    format     = "fusion_paged")
  arr <- array(as.double(seq_len(3 * 4 * 2)), dim = c(3, 4, 2),
               dimnames = list(NULL, NULL, c("DAPI", "CD3")))
  img <- bgnormR:::.new_QPTIFFImage(arr, meta)
  tmp <- tempfile(fileext = ".ome.tif")
  on.exit(unlink(tmp), add = TRUE)
  write_qptiff(img, tmp)

  m2 <- metadata(read_qptiff(tmp))
  expect_equal(qpi_pixel_size_um(m2), 0.5)
  expect_equal(m2$slide$slide_id, "SLIDE1")
  ct <- channel_table(m2)
  expect_equal(ct$fluorophore, c("DAPI", "OPAL520"))
  expect_equal(ct$exposure_time_us, c(3000, 150000))
  expect_equal(ct$color_rgb, c("0,0,255", "0,255,0"))
  # vendor-only field recovered from the qpi://vectra MapAnnotation
  expect_equal(qpi_channels(m2)[[1L]]$signal_units, 64L)
})

test_that("write_qptiff OME header declares channels (SizeC), not Z/T slices", {
  skip_if_not_installed("tiff")
  arr <- array(0, dim = c(3, 3, 4),
               dimnames = list(NULL, NULL, c("A", "B", "C", "D")))
  img <- as.QPTIFFImage(arr)
  tmp <- tempfile(fileext = ".ome.tif")
  on.exit(unlink(tmp), add = TRUE)
  write_qptiff(img, tmp)
  desc <- bgnormR:::.read_all_ifd_info(tmp)$descriptions[[1L]]
  expect_true(grepl('SizeC="4"', desc))
  expect_true(grepl('SizeZ="1"', desc))
  expect_true(grepl('SizeT="1"', desc))
})

# ---------- OME-Zarr reading (requires Rarr) ------------------------------

test_that(".match_channel_meta maps duplicate channel names in order", {
  meta <- bgnormR:::.single_scene_metadata(
    slide      = bgnormR:::.new_slide_info(),
    image_info = bgnormR:::.new_image_info(),
    channels   = list(list(index = 0L, name = "DAPI", fluorophore = "A"),
                      list(index = 1L, name = "CD3"),
                      list(index = 2L, name = "DAPI", fluorophore = "B")),
    scales     = bgnormR:::.build_scales(1L),
    format     = "fusion_paged")
  out <- bgnormR:::.match_channel_meta(meta, c("DAPI", "DAPI", "CD3"))
  expect_equal(out[[1L]]$fluorophore, "A")   # first DAPI -> first entry
  expect_equal(out[[2L]]$fluorophore, "B")   # second DAPI -> second entry
  expect_equal(out[[3L]]$name, "CD3")
})

test_that(".select_channels resolves names, indices, and errors on unknown", {
  nm <- c("DAPI", "CD3", "CD8")
  expect_equal(bgnormR:::.select_channels(NULL, nm)$ch_idx, 1:3)
  expect_equal(bgnormR:::.select_channels(c("CD8", "DAPI"), nm)$ch_idx, c(3L, 1L))
  expect_equal(bgnormR:::.select_channels(c(2L, 3L), nm)$ch_names, c("CD3", "CD8"))
  expect_error(bgnormR:::.select_channels("MISSING", nm), "not found")
})

test_that(".ome_page_indices falls back to identity for sparse TiffData", {
  ome <- list(
    sizes = list(size_c = 3L, size_z = 1L, size_t = 1L),
    dimension_order = "XYCZT",
    tiffdata = list(list(first_c = 0L, first_z = 0L, first_t = 0L,
                         ifd = 0L, plane_count = 1L)))  # only plane 0 mapped
  # channels 1,2 are absent from the table -> identity fallback, no crash
  expect_equal(bgnormR:::.ome_page_indices(ome, 3L), c(1L, 2L, 3L))
})

test_that(".is_ome_zarr_store detects .zarr paths and stores", {
  expect_true(bgnormR:::.is_ome_zarr_store("/some/where/img.ome.zarr"))
  expect_false(bgnormR:::.is_ome_zarr_store("/some/where/img.ome.tif"))
})

test_that("read_qptiff reads a synthetic OME-Zarr store", {
  skip_if_not_installed("Rarr")
  store <- file.path(tempdir(), paste0("zt_", as.integer(runif(1, 1, 1e6)), ".ome.zarr"))
  on.exit(unlink(store, recursive = TRUE), add = TRUE)
  dir.create(file.path(store, "scale0"), recursive = TRUE)
  a <- array(0L, dim = c(2, 5, 4))          # (c, y, x)
  a[1, , ] <- 10L; a[2, , ] <- matrix(1:20, 5, 4)
  Rarr::write_zarr_array(a, file.path(store, "scale0", "image"),
                         chunk_dim = c(1, 5, 4))
  Rarr::write_zarr_attributes(store, list(
    ome = list(
      version = "0.5",
      multiscales = list(list(
        axes = list(list(name = "c", type = "channel"),
                    list(name = "y", type = "space", unit = "micrometer"),
                    list(name = "x", type = "space", unit = "micrometer")),
        datasets = list(list(path = "scale0/image",
          coordinateTransformations = list(list(type = "scale",
                                                 scale = list(1, 0.5, 0.5))))),
        name = "syn")),
      omero = list(channels = list(list(label = "DAPI", color = "0000FF"),
                                   list(label = "CD3",  color = "FF0000"))))
  ))

  z <- read_qptiff(store)
  expect_equal(qpi_format(metadata(z)), "ome_zarr")
  expect_equal(names(z), c("DAPI", "CD3"))
  expect_equal(dim(z), c(5L, 4L, 2L))          # permuted (c,y,x) -> (y,x,c)
  expect_equal(qpi_pixel_size_um(metadata(z)), 0.5)
  expect_equal(as.numeric(z[1, 1, 1]), 10)     # channel 1
  expect_equal(as.numeric(z[5, 4, 2]), 20)     # channel 2, last pixel

  # lazy + channel subset
  zl <- read_qptiff(store, channels = "CD3", lazy = TRUE)
  expect_equal(names(zl), "CD3")
  expect_equal(as.numeric(as.array(zl[1, 1, 1, drop = FALSE])), 1)
})
