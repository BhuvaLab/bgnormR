# Classify observations using threshold-corrected GMM assignment

Starts from the posterior-probability argmax, then applies two
correction passes for the 3-component model so that the distribution
lower tails never overwhelm the hard evidence of a low pixel value:

1.  If argmax = 3 (Signal) but `x < means[2]` -\> demote to 2.

2.  If class is now 2 (Non-specific) but `x < means[1]` -\> demote to 1.

For the 2-component model the argmax is returned unchanged.

## Usage

``` r
.classify_pixels(x, means, posteriors)
```

## Arguments

- x:

  Numeric vector of (log-transformed) intensities, same length as
  `nrow(posteriors)`.

- means:

  Resolved component means (1 = background, 2 = non-specific, 3 =
  signal).

- posteriors:

  Matrix (n x G) of posterior probabilities.

## Value

Integer vector of component labels (values in `1:G`).
