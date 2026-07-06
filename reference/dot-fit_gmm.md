# Fit a Gaussian Mixture Model to intensity data

Wraps
[`mclust::Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html)
with the variable-variance (`"V"`) model and returns a tidy list of
parameters and posterior probabilities for the fitting subset. Use
[`.gmm_posteriors`](dot-gmm_posteriors.md) to obtain posteriors for the
full input vector.

## Usage

``` r
.fit_gmm(x, n_components = 3L, sample_prop = 1, .sample_threshold = 1e+05, ...)
```

## Arguments

- x:

  Numeric vector of (log-transformed) intensities.

- n_components:

  Integer, number of GMM components (2 or 3).

- sample_prop:

  Numeric in (0, 1\]; proportion of non-zero observations to use for
  fitting. Values \< 1 draw a random subset without replacement. Default
  `1` (use all non-zero observations).

- ...:

  Extra arguments forwarded to
  [`mclust::Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html).

## Value

A list with elements:

- means:

  Component means (length `n_components`).

- sds:

  Component standard deviations.

- props:

  Mixing proportions.

- posteriors:

  Matrix (n_fit x G) of posterior probabilities for the fitting subset
  only. Use [`.gmm_posteriors()`](dot-gmm_posteriors.md) for the full
  data.
