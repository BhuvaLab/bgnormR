# Background deconvolution for a GMM-fitted intensity vector

Implements the conditional-expectation deconvolution from Kharbanda et
al. 2025. Works for both the 3-component (pixel-level) and 2-component
(cell-level) models by using the last two components (G-1 and G) as the
background and signal pair:

## Usage

``` r
.deconvolve(x, means, sds, posteriors)
```

## Arguments

- x:

  Numeric vector of log-transformed intensities.

- means:

  Ordered GMM means (length G).

- sds:

  Ordered GMM standard deviations (length G).

- posteriors:

  Matrix (n x G) of posterior probabilities.

## Value

Numeric vector of background-adjusted log-intensities (same length as
`x`).

## Details

\$\$X\_{\text{adj}}(x) = P(C=G \mid x) \Bigl\[ (\mu_G - \mu\_{G-1}) +
\frac{\sigma_G^2 + \sigma\_{G-1}^2 - 2\min(\sigma\_{G-1}^2,\sigma_G^2)}
{\sigma_G^2} (x - \mu_G) \Bigr\]\$\$

For pixel-level (G=3): component 3 = signal, component 2 = non-specific.
For cell-level (G=2): component 2 = signal, component 1 = non-specific.
