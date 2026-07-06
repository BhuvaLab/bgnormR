# Read full page-layout metadata from every IFD in the TIFF chain

Makes one linear pass through all TIFF/BigTIFF IFDs and collects the
structural tags (width, height, compression, tile/strip offsets, etc.)
needed to read pixel data for each page later. No pixel data is read.

## Usage

``` r
.read_all_ifd_page_layouts(path, max_pages = 1000L)
```

## Arguments

- path:

  File path.

- max_pages:

  Upper bound on number of IFDs to traverse.

## Value

A list (one element per IFD) of named lists with fields: `image_width`,
`image_height`, `bits_per_sample`, `compression`, `tile_width`,
`tile_height`, `tile_offsets`, `tile_byte_counts`, `strip_offsets`,
`strip_byte_counts`, `rows_per_strip`, `sample_format`, `endian`.
