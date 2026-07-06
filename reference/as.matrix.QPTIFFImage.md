# Coerce a single-channel QPTIFFImage to a matrix

Drops the singleton channel dimension and returns a plain `H x W`
matrix. The channel name is stored in `attr(result, "channel")`. If the
object has more than one channel, use `img[, , "DAPI"]` to select one
first.

## Usage

``` r
# S3 method for class 'QPTIFFImage'
as.matrix(x, ...)
```

## Arguments

- x:

  A `QPTIFFImage` with exactly one channel, *or* a 2-D array (already a
  matrix-like object), which is returned with
  [`as.matrix`](https://rdrr.io/r/base/matrix.html) applied directly.

- ...:

  Unused.

## Value

A numeric matrix of dimension `c(H, W)`.

## Examples

``` r
arr <- array(1:12, dim = c(3, 4, 1),
             dimnames = list(NULL, NULL, "DAPI"))
img <- as.QPTIFFImage(arr)
m   <- as.matrix(img)
dim(m)          # 3 x 4
#> [1] 3 4
attr(m, "channel")  # "DAPI"
#> [1] "DAPI"
```
