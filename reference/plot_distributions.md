# Per-marker intensity distribution with fitted GMM densities

Plots a histogram of log2-transformed pixel (or cell) intensities and
overlays the component-wise density curves from the fitted GMM. One
facet per marker. Accepts a [`QPTIFFImage`](QPTIFFImage.md) returned by
[`bgnorm_pixels`](bgnorm_pixels.md) (bgnorm results are carried as an
attribute) or a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
/
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
with bgnorm results in `metadata(x)$bgnorm_results`. Multi-sample
`SpatialExperiment`s are automatically detected via `sample_id` in
`colData` and produce facets by `sample_id x marker`.

## Usage

``` r
plot_distributions(x, results = NULL, markers = NULL, ncol = NULL)
```

## Arguments

- x:

  A [`QPTIFFImage`](QPTIFFImage.md) returned by
  [`bgnorm_pixels`](bgnorm_pixels.md) *or* a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  /
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- results:

  Optionally, the [`QPTIFFImage`](QPTIFFImage.md) returned by
  [`bgnorm_pixels`](bgnorm_pixels.md) when `x` is the raw (unadjusted)
  `QPTIFFImage`. Ignored for `SummarizedExperiment` input.

- markers:

  Character vector of markers to display, or `NULL` for all markers.

- ncol:

  Number of facet columns.

## Value

A `ggplot` object.

## See also

[`plot_qptiff`](plot_qptiff.md),
[`plot_pixel_classes`](plot_pixel_classes.md),
[`bgnorm_pixels`](bgnorm_pixels.md), [`bgnorm_sce`](bgnorm_sce.md)

## Examples

``` r
path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
img  <- read_qptiff(path)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
res  <- bgnorm_pixels(img, sample_prop = 0.1)
plot_distributions(res)
```
