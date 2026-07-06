# Resolve GMM component identities

Orders components as background (1), non-specific/autofluorescence (2),
and biological signal (3). When the two highest-mean components are
close, the one with the heavier right tail (higher 75th percentile) is
designated signal, following Kharbanda et al. 2025.

## Usage

``` r
.resolve_components(means, sds, props, posteriors)
```

## Arguments

- means:

  Numeric vector of component means.

- sds:

  Numeric vector of component standard deviations.

- props:

  Numeric vector of mixing proportions.

- posteriors:

  Matrix of posterior probabilities (n x G).

## Value

A list with `means`, `sds`, `props`, and `posteriors` reordered so
component 1 = background, 2 = nonspecific, 3 = signal.
