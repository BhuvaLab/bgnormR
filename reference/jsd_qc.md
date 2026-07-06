# Compute the JSD quality control metric from bgnorm model parameters

Extracts the JSD between components 2 (non-specific binding) and 3
(signal) from the fitted model parameters stored in a
[`BgnormResult`](BgnormResult.md) object.

## Usage

``` r
jsd_qc(result)
```

## Arguments

- result:

  A `BgnormResult` object returned by
  [`bgnorm_pixels`](bgnorm_pixels.md) or
  [`bgnorm_cells`](bgnorm_cells.md).

## Value

Named numeric scalar `jsd`.

## Examples

``` r
set.seed(1)
x   <- c(exp(rnorm(300, log(300), 0.5)), exp(rnorm(100, log(2000), 0.9)))
res <- bgnorm_cells(x)
jsd_qc(res)
#> [1] 0.690225
```
