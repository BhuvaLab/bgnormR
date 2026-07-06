# Materialise a QPTIFFImage to a plain 3-D array

For eager QPTIFFImages, equivalent to `unclass(x)`. For lazy (on-disk)
QPTIFFImages, reads all pixel data from disk.

## Usage

``` r
# S3 method for class 'QPTIFFImage'
as.array(x, ...)
```

## Arguments

- x:

  A `QPTIFFImage`.

- ...:

  Unused.

## Value

A plain numeric array of dimension `c(H, W, C)`.
