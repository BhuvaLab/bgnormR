# Spatial GMM class assignment plot for a QPTIFFImage or SpatialExperiment

Assigns each pixel/cell to its most probable GMM component and renders
the result. For a [`QPTIFFImage`](QPTIFFImage.md) produced by
[`bgnorm_pixels`](bgnorm_pixels.md): pseudocoloured raster (pixel-level
data, three components). For a
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html):
scatter plot at spatial coordinates using the bgnorm results from
`metadata(spe)$bgnorm_results` (two components). Multi-sample
`SpatialExperiment`s are automatically detected via the `sample_id`
column in `colData` and faceted accordingly.

## Usage

``` r
plot_pixel_classes(
  x,
  markers = NULL,
  resolution = 1L,
  point_size = 1,
  pixels = c(1024L, 1024L),
  flip_y = TRUE,
  large_data_threshold = 10000L,
  ncol = NULL
)
```

## Arguments

- x:

  A [`QPTIFFImage`](QPTIFFImage.md) returned by
  [`bgnorm_pixels`](bgnorm_pixels.md) (contains both the adjusted pixel
  intensities and the embedded [`BgnormResult`](BgnormResult.md)
  objects), *or* a
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  with `metadata(x)$bgnorm_results` populated by
  [`bgnorm_sce`](bgnorm_sce.md).

- markers:

  Character vector of markers to display, or `NULL` for all markers in
  `x`.

- resolution:

  Positive integer; pixel downsample factor for raster output. Ignored
  for `SpatialExperiment`.

- point_size:

  Base point size for `SpatialExperiment` scatter.

- pixels:

  Integer `c(width, height)` for scattermore rasterisation. Default
  `c(1024L, 1024L)`.

- flip_y:

  Logical; reverse y-axis for `SpatialExperiment`? Default `TRUE`.

- large_data_threshold:

  Cell count threshold for automatic scattermore use. Default `10000L`.

- ncol:

  Number of facet columns.

## Value

A `ggplot` object.

## See also

[`plot_qptiff`](plot_qptiff.md),
[`plot_distributions`](plot_distributions.md),
[`bgnorm_pixels`](bgnorm_pixels.md), [`bgnorm_sce`](bgnorm_sce.md)

## Examples

``` r
path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
img  <- read_qptiff(path)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
res  <- bgnorm_pixels(img, sample_prop = 0.1)
plot_pixel_classes(res, markers = c("PanCK", "CD20"))
```
