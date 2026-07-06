# Apply a 3 x 3 median filter to a matrix

Uses
[`EBImage::medianFilter`](https://rdrr.io/pkg/EBImage/man/medianFilter.html)
(C-based, fast) when EBImage is available, otherwise falls back to a
pure-R sliding-window implementation.

## Usage

``` r
median_filter_3x3(img)
```

## Arguments

- img:

  Numeric matrix (height x width).

## Value

Filtered matrix of the same dimensions.

## Examples

``` r
m <- matrix(runif(100), 10, 10)
m_filt <- median_filter_3x3(m)
```
