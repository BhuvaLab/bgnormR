# Inverse log2 transform with cofactor

Reverses [`log_transform`](log_transform.md): \\I = (2^{I\_{log}} - 1)
\times \text{cofactor}\\.

## Usage

``` r
inv_log_transform(x, cofactor = 150)
```

## Arguments

- x:

  Numeric vector or matrix of log2-transformed values.

- cofactor:

  Positive numeric cofactor used during `log_transform`.

## Value

Numeric vector or matrix on the original intensity scale.

## Examples

``` r
x <- c(0, 100, 500, 2000, 65535)
x_log <- log_transform(x)
all.equal(x, inv_log_transform(x_log))
#> [1] TRUE
```
