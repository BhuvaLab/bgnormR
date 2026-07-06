# Decode a single compressed tile or strip via an in-memory TIFF wrapper

Wraps the raw compressed bytes in a minimal single-strip TIFF structure
and decodes via
[`tiff::readTIFF()`](https://rdrr.io/pkg/tiff/man/readTIFF.html). This
delegates all decompression (LZW, JPEG, PackBits, etc.) to libtiff
without loading the entire source file.

## Usage

``` r
.decode_tile_via_tiff(raw_bytes, w, h, bps, compression, endian = "little")
```

## Arguments

- raw_bytes:

  Raw vector of compressed tile/strip bytes.

- w, h:

  Tile/strip width and height in pixels.

- bps:

  Bits per sample (8, 16, or 32).

- compression:

  TIFF compression code.

- endian:

  Byte order for the in-memory TIFF header.

## Value

Integer matrix `[h, w]`, or `NULL` on error.
