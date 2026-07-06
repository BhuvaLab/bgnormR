# Read tag 270 and tag 277 from every IFD in a TIFF/BigTIFF file

Read tag 270 and tag 277 from every IFD in a TIFF/BigTIFF file

## Usage

``` r
.read_all_ifd_info(path, max_pages = 1000L)
```

## Arguments

- path:

  File path.

- max_pages:

  Maximum IFDs to traverse (default 1000).

## Value

List with `descriptions` (character) and `samples_per_pixel` (integer).
