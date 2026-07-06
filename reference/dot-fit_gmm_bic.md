# Fit G=2 and G=3 GMMs on the same pixel sample and compare BIC

Filters zeros, optionally subsamples, then calls
[`mclust::Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html)
with `G = 2:3` so that both models are fitted on an identical dataset
and their BIC values are directly comparable. mclust selects the
higher-BIC model; `no_signal` is `TRUE` when G=2 is selected.

## Usage

``` r
.fit_gmm_bic(x, sample_prop = 1, .sample_threshold = 1e+05, ...)
```

## Arguments

- x:

  Numeric vector of log-transformed intensities.

- sample_prop:

  Numeric in (0, 1\]; subsampling fraction (default 1).

- .sample_threshold:

  Integer; only subsample when `n_nonzero` exceeds this value (default
  1e5).

- ...:

  Extra arguments forwarded to
  [`mclust::Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html).

## Value

A list with elements `means`, `sds`, `props`, `posteriors` (fitting
subset, G components), `bic` (named `c(G2 = ..., G3 = ...)`, each
divided by the number of points used to fit so values are point-specific
and comparable across channels), and `no_signal` (logical).
