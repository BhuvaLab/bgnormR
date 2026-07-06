# Compute the model-based quantile normalisation factor

Evaluates the `quantile`th quantile of the signal component distribution
and maps it through the bgnorm deconvolution formula to obtain a
normalisation factor.

## Usage

``` r
.qnorm_factor(means, sds, props, quantile = 0.75)
```

## Arguments

- means:

  Ordered component means (length G; G=3 for pixel, G=2 for cell).

- sds:

  Ordered component standard deviations.

- props:

  Mixing proportions (from the fitted GMM).

- quantile:

  Quantile of the signal component to use (default 0.75).

## Value

Scalar normalisation factor (positive).
