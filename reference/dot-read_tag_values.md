# Read a TIFF tag's value(s) as a numeric vector

Handles inline vs offset storage and all common TIFF types (BYTE, SHORT,
LONG, LONG8). Returns a double vector for large offsets to avoid signed
32-bit integer overflow.

## Usage

``` r
.read_tag_values(f, tiff_type, count, vfield, bigtiff, endian)
```
