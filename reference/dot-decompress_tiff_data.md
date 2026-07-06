# Decompress one TIFF tile or strip

Handles uncompressed (compression = 1) and DEFLATE/zlib (8, 32946). For
DEFLATE, uses `memDecompress(type="gzip")` which correctly handles
zlib-wrapped deflate streams (0x78 0x9C / 0x78 0xDA headers) produced by
the Akoya PhenoCycler-Fusion software.

## Usage

``` r
.decompress_tiff_data(raw_bytes, compression, n_bytes_out)
```

## Arguments

- raw_bytes:

  Raw vector of compressed bytes from file.

- compression:

  TIFF compression code.

- n_bytes_out:

  Expected number of bytes after decompression.

## Value

Raw vector of decompressed bytes.
