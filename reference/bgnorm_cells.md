# Cell-level background normalisation (bgnorm)

Applies a two-component Gaussian Mixture Model to cell-level aggregated
intensities for a single marker. The two components represent
non-specific binding / autofluorescence and biological signal.

## Usage

``` r
bgnorm_cells(
  x,
  cofactor = 150,
  quantile_norm = FALSE,
  quantile = 0.75,
  sample_prop = 1,
  ...
)
```

## Arguments

- x:

  Numeric vector of raw cell-level intensities for one marker. Must be a
  plain vector; matrices and arrays are not accepted.

- cofactor:

  Positive numeric cofactor for log2 transform (default 150).

- quantile_norm:

  Logical; apply model-based quantile normalisation?

- quantile:

  Quantile of the signal component for normalisation (default 0.75).

- sample_prop:

  Numeric in (0, 1\]; proportion of non-zero observations used to fit
  the GMM. Default `1` (use all).

- ...:

  Additional arguments forwarded to the internal GMM fitter.

## Value

A `BgnormResult` object.

## Details

This is a cell-level approximation of the pixel-level method intended
for cases where only cell-summarised intensities (e.g., mean or median
intensity per cell) are available.

## See also

[`bgnorm_pixels`](bgnorm_pixels.md), [`bgnorm_sce`](bgnorm_sce.md)

## Examples

``` r
set.seed(3)
# Simulate cell intensities: nonspecific + signal
x <- c(exp(rnorm(300, log(200), 0.5)), exp(rnorm(100, log(1500), 0.8)))
res <- bgnorm_cells(x)
print(res)
#> BgnormResult (cell-level)
#>   n = 400 
#>   Component means: 1.236 3.139 
#>   JSD (QC metric): 0.6575 
#>   No signal detected: FALSE 
#>   Quantile normalised: FALSE 
#>   Tissue positivity: 31.3% 
```
