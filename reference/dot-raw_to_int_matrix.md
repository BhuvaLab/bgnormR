# Convert raw pixel bytes to an integer matrix

Convert raw pixel bytes to an integer matrix

## Usage

``` r
.raw_to_int_matrix(raw_bytes, nr, nc, bps, endian)
```

## Arguments

- raw_bytes:

  Raw vector (decompressed pixel data, row-major).

- nr:

  Number of rows.

- nc:

  Number of columns.

- bps:

  Bits per sample (8 or 16).

- endian:

  Byte order ("little" or "big").

## Value

Integer matrix `[nr, nc]`.
