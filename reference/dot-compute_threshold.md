# Compute the Non-specific / Signal class boundary threshold

Given a GMM-based per-pixel classification (from `.classify_pixels`) and
the corresponding background-adjusted intensities, computes the maximum
adjusted value among pixels classified as Non-specific (class 2). This
single scalar is sufficient to recover the full 3-class assignment:

- Class 1 (Background): `adjusted == 0`

- Class 2 (Non-specific): `0 < adjusted <= threshold`

- Class 3 (Signal): `adjusted > threshold`

The class-1 boundary is always 0 (enforced by the `pmax` in
`.deconvolve`) and is not stored.

## Usage

``` r
.compute_threshold(x_adj, cls)
```

## Arguments

- x_adj:

  Numeric vector of background-adjusted intensities (length n, all
  values `>= 0` after `pmax`).

- cls:

  Integer vector of class labels produced by `.classify_pixels` (values
  in `1:3`).

## Value

A single numeric scalar: the maximum adjusted intensity of class-2
pixels, or `0` if no pixels are classified as Non-specific (in which
case every non-zero pixel is Signal by convention).
