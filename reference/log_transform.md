# Log2 transform with cofactor

Applies \\I\_{log} = \log_2(I / \text{cofactor} + 1)\\ to stabilise
variance in low-intensity fluorescence signals. The default cofactor of
150 is appropriate for 16-bit images from the Akoya PhenoCycler-Fusion;
adjust for other platforms.

## Usage

``` r
log_transform(x, cofactor = 150)
```

## Arguments

- x:

  Numeric vector or matrix of raw intensities.

- cofactor:

  Positive numeric cofactor (default 150 for 16-bit images).

## Value

Numeric vector or matrix of log2-transformed values.

## See also

[`inv_log_transform`](inv_log_transform.md)

## Examples

``` r
x <- c(0, 100, 500, 2000, 10000, 65535)
log_transform(x)
#> [1] 0.0000000 0.7369656 2.1154772 3.8413023 6.0803734 8.7744576
```
