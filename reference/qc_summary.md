# Summarise QC across multiple bgnorm results

Summarise QC across multiple bgnorm results

## Usage

``` r
qc_summary(results)
```

## Arguments

- results:

  A named list of `BgnormResult` objects (e.g., one per marker), or a
  `QPTIFFImage` returned by [`bgnorm_pixels`](bgnorm_pixels.md) (bgnorm
  results are extracted automatically).

## Value

A `data.frame` with columns `name`, `jsd`, and `prop_signal` (proportion
of pixels/cells in the signal component).

## Examples

``` r
set.seed(2)
mk1 <- bgnorm_cells(c(exp(rnorm(300,log(300),0.5)), exp(rnorm(100,log(2000),0.9))))
mk2 <- bgnorm_cells(c(exp(rnorm(400,log(250),0.4)), exp(rnorm( 80,log(1800),0.8))))
qc_summary(list(CD3 = mk1, CD8 = mk2))
#>   name       jsd prop_signal
#> 1  CD3 0.6863340   0.3158777
#> 2  CD8 0.7329585   0.1790243
```
