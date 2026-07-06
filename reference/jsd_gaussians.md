# Jensen-Shannon Divergence between two Gaussian distributions

Computes JSD numerically on a dense grid, used as a signal-to-noise QC
metric between the non-specific binding component (2) and the biological
signal component (3) of the bgnorm model. Higher values indicate better
separation and staining quality.

## Usage

``` r
jsd_gaussians(mu1, sd1, mu2, sd2, n_grid = 2000L)
```

## Arguments

- mu1, mu2:

  Means of the two Gaussian distributions.

- sd1, sd2:

  Standard deviations of the two Gaussian distributions.

- n_grid:

  Number of evaluation points (default 2000).

## Value

Scalar JSD value in \[0, 1\] (bits, log2 base).

## Examples

``` r
# Well-separated distributions -> high JSD
jsd_gaussians(mu1 = 0, sd1 = 0.5, mu2 = 3, sd2 = 0.8)
#> [1] 0.9577609
# Identical distributions -> JSD = 0
jsd_gaussians(mu1 = 1, sd1 = 1, mu2 = 1, sd2 = 1)
#> [1] 0
```
