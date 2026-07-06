# Compute GMM posterior probabilities for an arbitrary input vector

Given fitted GMM parameters (from [`.fit_gmm`](dot-fit_gmm.md)),
computes the posterior probability of each observation belonging to each
component. This should be called on the full data vector after fitting
on a subsample.

## Usage

``` r
.gmm_posteriors(x, means, props, sds)
```

## Arguments

- x:

  Numeric vector; the full set of observations.

- means:

  Component means.

- props:

  Mixing proportions.

- sds:

  Component standard deviations.

## Value

A matrix of dimension `length(x) x G` with rows summing to 1.
