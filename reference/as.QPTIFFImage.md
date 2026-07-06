# Coerce a 3-D array to a QPTIFFImage

Attaches the `QPTIFFImage` class to any 3-D numeric array `[H x W x C]`.
Channel names are taken from `dimnames(x)[[3]]`. This is the intended
way to promote a plain array produced outside `read_qptiff` (e.g., the
adjusted output of `bgnorm_markers`) into a first-class `QPTIFFImage`.

## Usage

``` r
as.QPTIFFImage(x, ...)

# S3 method for class 'QPTIFFImage'
as.QPTIFFImage(x, ...)

# S3 method for class 'array'
as.QPTIFFImage(x, metadata = list(), ...)

# S3 method for class 'DelayedArray'
as.QPTIFFImage(x, metadata = list(), ...)
```

## Arguments

- x:

  A 3-D numeric array `[H x W x C]`. For a `QPTIFFImage`, the object is
  returned unchanged.

- ...:

  Unused; present for S3 method consistency.

- metadata:

  Optional named list of metadata to embed. Ignored when `x` is already
  a `QPTIFFImage`.

## Value

A `QPTIFFImage`.

## Examples

``` r
arr <- array(runif(20 * 20 * 3), dim = c(20, 20, 3),
             dimnames = list(NULL, NULL, c("DAPI", "CD3", "CD8")))
img <- as.QPTIFFImage(arr)
class(img)
#> [1] "QPTIFFImage" "array"      
```
