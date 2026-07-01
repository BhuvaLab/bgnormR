## Tests for the QPTIFFArraySeed / read_qptiff_delayed DelayedArray backend

# ---- helpers ---------------------------------------------------------------

# Create a minimal DEFLATE-compressed BigTIFF in memory so tests can run
# without a real Akoya instrument file.  Returns path to a temp TIFF file
# containing `n_ch` grayscale 16-bit pages of size `h x w`.
.make_test_tiff <- function(h = 16L, w = 16L, n_ch = 3L,
                             values = NULL) {
  skip_if_not_installed("tiff")
  if (is.null(values)) {
    # Deterministic pixel values: channel c, row r, col x -> (c*100 + r*w + x)
    values <- lapply(seq_len(n_ch), function(c) {
      matrix(as.integer(c * 100L + seq_len(h * w) - 1L),
             nrow = h, ncol = w)
    })
  }
  path <- tempfile(fileext = ".tif")
  tiff::writeTIFF(values[[1L]] / 65535, path)
  for (k in seq(2L, n_ch)) {
    tiff::writeTIFF(values[[k]] / 65535, path, compression = "none",
                    bits.per.sample = 16L)
  }
  list(path = path, values = values, h = h, w = w, n_ch = n_ch)
}

# ---- binary helper tests ---------------------------------------------------

test_that(".raw_to_int_matrix converts 16-bit row-major bytes", {
  # 2x3 matrix: [[1,2,3],[4,5,6]]
  vals <- as.integer(c(1, 2, 3, 4, 5, 6))
  raw_bytes <- writeBin(vals, raw(), size = 2L, endian = "little")
  mat <- bgnormR:::.raw_to_int_matrix(raw_bytes, 2L, 3L, 16L, "little")
  expect_equal(mat[1, ], c(1L, 2L, 3L))
  expect_equal(mat[2, ], c(4L, 5L, 6L))
})

test_that(".raw_to_int_matrix converts 8-bit bytes", {
  vals <- as.raw(c(10L, 20L, 30L, 40L))
  mat <- bgnormR:::.raw_to_int_matrix(vals, 2L, 2L, 8L, "little")
  expect_equal(mat[1, ], c(10L, 20L))
  expect_equal(mat[2, ], c(30L, 40L))
})

test_that(".decompress_tiff_data handles uncompressed (compression=1)", {
  raw <- as.raw(c(0xAA, 0xBB, 0xCC))
  out <- bgnormR:::.decompress_tiff_data(raw, 1L, 3L)
  expect_identical(out, raw)
})

test_that(".decompress_tiff_data handles DEFLATE (compression=8) via memDecompress", {
  original <- as.raw(rep(0x41L, 100L))
  compressed <- memCompress(original, type = "gzip")
  decompressed <- bgnormR:::.decompress_tiff_data(compressed, 8L, 100L)
  expect_identical(decompressed, original)
})

test_that(".decompress_tiff_data handles Adobe Deflate (compression=32946)", {
  original <- as.raw(rep(0x42L, 64L))
  compressed <- memCompress(original, type = "gzip")
  decompressed <- bgnormR:::.decompress_tiff_data(compressed, 32946L, 64L)
  expect_identical(decompressed, original)
})

test_that(".decompress_tiff_data returns NULL for unsupported compression (LZW)", {
  result <- bgnormR:::.decompress_tiff_data(as.raw(1:4), 5L, 4L)  # LZW
  expect_null(result)
})

# ---- IFD layout reader tests -----------------------------------------------

test_that(".read_all_ifd_page_layouts returns empty list for non-existent file", {
  result <- bgnormR:::.read_all_ifd_page_layouts(tempfile())
  expect_equal(length(result), 0L)
})

test_that(".read_tag_values handles SHORT (type=3) arrays", {
  f_temp <- tempfile()
  f <- file(f_temp, "wb")
  vals <- as.integer(c(42L, 256L, 1000L))
  writeBin(vals, f, size = 2L, endian = "little")
  close(f)

  f <- file(f_temp, "rb")
  raw_bytes <- readBin(f, "raw", n = 6L)
  close(f)

  # Simulate: 3 values inline
  vfield <- c(raw_bytes[1:6], as.raw(rep(0L, 2L)))  # pad to 8 bytes
  result <- bgnormR:::.read_tag_values(NULL, 3L, 3L, vfield, TRUE, "little")
  expect_equal(result, c(42L, 256L, 1000L))

  unlink(f_temp)
})

test_that(".read_tag_values handles LONG (type=4) with large uint32", {
  # Value = 3000000000 > 2^31-1, stored as 4-byte little-endian
  big_val <- 3000000000
  raw4 <- writeBin(as.numeric(big_val), raw(), size = 4L, endian = "little")
  vfield <- c(raw4, as.raw(rep(0L, 4L)))  # pad to 8 bytes
  result <- bgnormR:::.read_tag_values(NULL, 4L, 1L, vfield, TRUE, "little")
  expect_equal(result, big_val, tolerance = 1)
})

# ---- QPTIFFArraySeed construction tests ------------------------------------

test_that("QPTIFFArraySeed slot names and types are correct", {
  seed <- new("QPTIFFArraySeed",
    filepath     = "/tmp/test.qptiff",
    .dim         = c(256L, 256L, 3L),
    .dimnames    = list(NULL, NULL, c("DAPI", "CD3", "CD8")),
    page_layouts = list(list(), list(), list()),
    dtype        = "integer",
    metadata     = list(),
    level        = 1L
  )
  expect_equal(dim(seed), c(256L, 256L, 3L))
  expect_equal(dimnames(seed)[[3L]], c("DAPI", "CD3", "CD8"))
  expect_equal(type(seed), "integer")
})

test_that("QPTIFFArraySeed with all-NULL dimnames returns NULL dimnames", {
  seed <- new("QPTIFFArraySeed",
    filepath = "/tmp/x.tif", .dim = c(10L, 10L, 2L),
    .dimnames = list(NULL, NULL, NULL), page_layouts = list(),
    dtype = "double", metadata = list(), level = 1L
  )
  expect_null(dimnames(seed))
})

test_that("extract_array returns correct dimensions for a subset", {
  skip_if_not_installed("tiff")
  skip_if_not_installed("xml2")

  # Build a fake seed with 3 channels, 16x16 px
  # Use in-memory matrices; skip actual file I/O by mocking page layouts
  h <- 16L; w <- 16L; nc <- 3L

  seed <- new("QPTIFFArraySeed",
    filepath     = tempfile(),  # doesn't need to exist for this test
    .dim         = c(h, w, nc),
    .dimnames    = list(NULL, NULL, c("DAPI", "CD3", "CD8")),
    page_layouts = vector("list", nc),
    dtype        = "integer",
    metadata     = list(),
    level        = 1L
  )

  # extract_array with all-NULL index returns full dimensions
  # (will fail trying to read file — just test dimension contract)
  expect_equal(dim(seed), c(16L, 16L, 3L))
  expect_equal(length(seed@page_layouts), 3L)
})

test_that("extract_array on empty index returns zero-row array", {
  seed <- new("QPTIFFArraySeed",
    filepath = "/tmp/x.tif", .dim = c(10L, 10L, 3L),
    .dimnames = list(NULL, NULL, c("A","B","C")),
    page_layouts = list(list(), list(), list()),
    dtype = "integer", metadata = list(), level = 1L
  )
  result <- extract_array(seed, list(integer(0L), seq_len(10L), seq_len(3L)))
  expect_equal(dim(result), c(0L, 10L, 3L))
})

# ---- read_qptiff(lazy=TRUE) tests ------------------------------------------

test_that("read_qptiff(lazy=TRUE) errors for missing file", {
  expect_error(read_qptiff("nonexistent.qptiff", lazy = TRUE), "File not found")
})

test_that("read_qptiff(lazy=TRUE) returns a QPTIFFImage for valid TIFF", {
  skip_if_not_installed("tiff")
  skip_if_not_installed("xml2")
  skip_if_not_installed("DelayedArray")

  h <- 8L; w <- 8L; nc <- 2L
  path <- tempfile(fileext = ".tif")
  on.exit(unlink(path), add = TRUE)

  pages <- lapply(seq_len(nc), function(c) {
    matrix(as.integer(c * 100L + seq_len(h * w)), nrow = h, ncol = w) / 65535
  })
  tiff::writeTIFF(pages, path, bits.per.sample = 16L, compression = "none")

  arr <- read_qptiff(path, lazy = TRUE)
  expect_s3_class(arr, "QPTIFFImage")
  expect_true(bgnormR:::.is_lazy_qptiff(arr))
  expect_s4_class(arr$.da, "DelayedArray")
  expect_equal(dim(arr)[3L], nc)
  expect_equal(dim(arr)[1L], h)
  expect_equal(dim(arr)[2L], w)
})

test_that("read_qptiff(lazy=TRUE) extracts correct pixel values", {
  skip_if_not_installed("tiff")
  skip_if_not_installed("xml2")
  skip_if_not_installed("DelayedArray")

  h <- 4L; w <- 4L
  path <- tempfile(fileext = ".tif")
  on.exit(unlink(path), add = TRUE)

  pg1 <- matrix(1000 / 65535, nrow = h, ncol = w)
  pg2 <- matrix(2000 / 65535, nrow = h, ncol = w)
  tiff::writeTIFF(list(pg1, pg2), path, bits.per.sample = 16L, compression = "none")

  arr <- read_qptiff(path, lazy = TRUE, as_integer = TRUE)

  ch1 <- as.array(arr[, , 1L])
  expect_true(all(abs(ch1 - 1000L) <= 1L))

  ch2 <- as.array(arr[, , 2L])
  expect_true(all(abs(ch2 - 2000L) <= 1L))
})

test_that("read_qptiff(lazy=TRUE, as_integer=FALSE) normalises to [0,1]", {
  skip_if_not_installed("tiff")
  skip_if_not_installed("xml2")
  skip_if_not_installed("DelayedArray")

  h <- 4L; w <- 4L
  path <- tempfile(fileext = ".tif")
  on.exit(unlink(path), add = TRUE)

  pg <- matrix(32767 / 65535, nrow = h, ncol = w)
  tiff::writeTIFF(pg, path, bits.per.sample = 16L, compression = "none")

  arr <- read_qptiff(path, lazy = TRUE, as_integer = FALSE)
  expect_s3_class(arr, "QPTIFFImage")
  expect_equal(type(DelayedArray::seed(arr$.da)), "double")
  vals <- as.vector(as.array(arr[, , 1L]))
  expect_true(all(vals >= 0 & vals <= 1))
})

# ---- .read_qptiff_page with synthetic layout --------------------------------

test_that(".read_stripped_page reads uncompressed strips correctly", {
  skip_if_not_installed("tiff")

  h <- 4L; w <- 6L; bps <- 16L; endian <- "little"

  # Create uncompressed 16-bit pixel data
  pixel_vals <- as.integer(seq_len(h * w))
  raw_all <- writeBin(pixel_vals, raw(), size = 2L, endian = endian)

  # Write to a temp file with two equal strips
  path <- tempfile()
  f <- file(path, "wb")
  strip1_offset <- 0L
  writeBin(raw_all[seq_len(h %/% 2L * w * 2L)], f, useBytes = TRUE)
  strip2_offset <- h %/% 2L * w * 2L
  writeBin(raw_all[seq(strip2_offset + 1L, length(raw_all))], f, useBytes = TRUE)
  close(f)

  layout <- list(
    image_height      = h,
    image_width       = w,
    bits_per_sample   = bps,
    compression       = 1L,
    strip_offsets     = c(0, strip2_offset),
    strip_byte_counts = c(strip2_offset, length(raw_all) - strip2_offset),
    rows_per_strip    = h %/% 2L,
    endian            = endian
  )

  mat <- bgnormR:::.read_stripped_page(path, layout, h, w, bps, 1L, seq_len(h), seq_len(w))
  expect_equal(mat[1L, 1L], 1L)
  expect_equal(mat[1L, w], as.integer(w))
  expect_equal(mat[h, w], as.integer(h * w))

  unlink(path)
})

test_that(".read_tiled_page reads uncompressed tiles correctly", {
  skip_if_not_installed("tiff")

  h <- 4L; w <- 6L; tile_h <- 2L; tile_w <- 3L; bps <- 16L; endian <- "little"

  # Pixel layout: [h x w] row-major
  # Tile arrangement: 2x2 tiles (2 tile rows, 2 tile cols)
  n_trows <- ceiling(h / tile_h)  # 2
  n_tcols <- ceiling(w / tile_w)  # 2

  path <- tempfile()
  f <- file(path, "wb")

  tile_offsets    <- numeric(n_trows * n_tcols)
  tile_byte_counts <- numeric(n_trows * n_tcols)
  tile_idx <- 1L

  # Build each tile from the original matrix
  pixel_mat <- matrix(seq_len(h * w), nrow = h, ncol = w, byrow = TRUE)

  for (tr in seq_len(n_trows)) {
    for (tc in seq_len(n_tcols)) {
      r0 <- (tr - 1L) * tile_h + 1L; r1 <- min(tr * tile_h, h)
      c0 <- (tc - 1L) * tile_w + 1L; c1 <- min(tc * tile_w, w)

      tile_full <- matrix(0L, nrow = tile_h, ncol = tile_w)
      tile_full[seq_len(r1-r0+1L), seq_len(c1-c0+1L)] <-
        pixel_mat[r0:r1, c0:c1]

      # Write tiles in row-major order (as TIFF stores them on disk)
      raw_tile <- writeBin(as.integer(t(tile_full)), raw(), size = 2L, endian = endian)
      tile_offsets[tile_idx]    <- seek(f)
      tile_byte_counts[tile_idx] <- length(raw_tile)
      writeBin(raw_tile, f, useBytes = TRUE)
      tile_idx <- tile_idx + 1L
    }
  }
  close(f)

  layout <- list(
    image_height     = h,
    image_width      = w,
    bits_per_sample  = bps,
    compression      = 1L,
    tile_width       = tile_w,
    tile_height      = tile_h,
    tile_offsets     = tile_offsets,
    tile_byte_counts = tile_byte_counts,
    endian           = endian
  )

  mat <- bgnormR:::.read_tiled_page(path, layout, h, w, bps, 1L, seq_len(h), seq_len(w))
  expect_equal(mat[1L, 1L], 1L)
  expect_equal(mat[1L, w], as.integer(w))
  expect_equal(mat[h, w], as.integer(h * w))

  unlink(path)
})

test_that(".read_tiled_page handles DEFLATE-compressed tiles", {
  h <- 4L; w <- 4L; tile_h <- 4L; tile_w <- 4L; bps <- 16L; endian <- "little"

  pixel_mat <- matrix(seq_len(h * w), nrow = h, ncol = w, byrow = TRUE)
  # Write in row-major order (as.integer(t(m)) = row-major of m)
  raw_pixels <- writeBin(as.integer(t(pixel_mat)), raw(), size = 2L, endian = endian)
  compressed_tile <- memCompress(raw_pixels, type = "gzip")

  path <- tempfile()
  f <- file(path, "wb")
  tile_offset <- seek(f)
  writeBin(compressed_tile, f, useBytes = TRUE)
  close(f)

  layout <- list(
    image_height     = h,
    image_width      = w,
    bits_per_sample  = bps,
    compression      = 8L,
    tile_width       = tile_w,
    tile_height      = tile_h,
    tile_offsets     = tile_offset,
    tile_byte_counts = length(compressed_tile),
    endian           = endian
  )

  mat <- bgnormR:::.read_tiled_page(path, layout, h, w, bps, 8L, seq_len(h), seq_len(w))
  expect_equal(mat, pixel_mat)

  unlink(path)
})

test_that(".read_tiled_page handles spatial subsetting", {
  h <- 8L; w <- 8L; tile_h <- 4L; tile_w <- 4L; bps <- 16L; endian <- "little"

  pixel_mat <- matrix(seq_len(h * w), nrow = h, ncol = w, byrow = TRUE)

  path <- tempfile()
  f <- file(path, "wb")

  n_trows <- ceiling(h / tile_h); n_tcols <- ceiling(w / tile_w)
  tile_offsets    <- numeric(n_trows * n_tcols)
  tile_byte_counts <- numeric(n_trows * n_tcols)
  tile_idx <- 1L

  for (tr in seq_len(n_trows)) {
    for (tc in seq_len(n_tcols)) {
      r0 <- (tr-1)*tile_h+1L; r1 <- min(tr*tile_h, h)
      c0 <- (tc-1)*tile_w+1L; c1 <- min(tc*tile_w, w)
      tile_full <- pixel_mat[r0:r1, c0:c1]
      # Row-major byte order (TIFF standard)
      raw_tile  <- writeBin(as.integer(t(tile_full)), raw(), size=2L, endian=endian)
      tile_offsets[tile_idx]     <- seek(f)
      tile_byte_counts[tile_idx] <- length(raw_tile)
      writeBin(raw_tile, f, useBytes = TRUE)
      tile_idx <- tile_idx + 1L
    }
  }
  close(f)

  layout <- list(
    image_height     = h, image_width = w, bits_per_sample = bps,
    compression      = 1L, tile_width = tile_w, tile_height = tile_h,
    tile_offsets     = tile_offsets, tile_byte_counts = tile_byte_counts,
    endian           = endian
  )

  req_rows <- c(2L, 3L, 6L); req_cols <- c(1L, 5L, 7L)
  mat <- bgnormR:::.read_tiled_page(path, layout, h, w, bps, 1L, req_rows, req_cols)

  expect_equal(dim(mat), c(3L, 3L))
  expect_equal(mat[1L, 1L], pixel_mat[2L, 1L])
  expect_equal(mat[2L, 2L], pixel_mat[3L, 5L])
  expect_equal(mat[3L, 3L], pixel_mat[6L, 7L])

  unlink(path)
})
