# Precompute a histogram of log-transformed intensities for storage

Used during `bgnorm_pixels` and `bgnorm_cells` fitting to store a
compact histogram of the raw data so that `plot_distributions` does not
need the original image or assay.

## Usage

``` r
.compute_histogram(x_log, bins = 100L)
```

## Arguments

- x_log:

  Numeric vector of log-transformed intensities (may contain zeros and
  non-finite values; these are excluded before binning).

- bins:

  Integer number of histogram bins (default 100).

## Value

A list with elements `breaks` (length bins + 1) and `density` (length
bins), as returned by `hist(..., plot = FALSE)`.
