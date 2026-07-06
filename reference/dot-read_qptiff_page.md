# Read one TIFF page as a numeric matrix

Routes to the appropriate reader based on the page layout: tiled
DEFLATE, tiled uncompressed, stripped DEFLATE, stripped uncompressed, or
a fallback through
[`tiff::readTIFF`](https://rdrr.io/pkg/tiff/man/readTIFF.html).

## Usage

``` r
.read_qptiff_page(path, layout, rows = NULL, cols = NULL)
```

## Arguments

- path:

  File path.

- layout:

  Page layout list from `.read_all_ifd_page_layouts`.

- rows:

  Integer vector of 1-based row indices to return, or `NULL` for all
  rows.

- cols:

  Integer vector of 1-based column indices to return, or `NULL` for all
  columns.

## Value

Integer matrix `[length(rows), length(cols)]`.
