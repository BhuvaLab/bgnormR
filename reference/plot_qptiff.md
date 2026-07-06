# Spatial intensity plot for a QPTIFFImage or SpatialExperiment

For a [`QPTIFFImage`](QPTIFFImage.md): renders up to 10 channels as an
additive colour composite on a black background. Each channel is
assigned a distinct colour; its per-pixel intensity (min-max normalised
to \[0, 1\]) drives the channel's contribution - high-intensity pixels
appear fully saturated while low-intensity pixels are transparent
(black). Channels that overlap in space produce mixed additive colours
(e.g. cyan + red -\> white), matching the standard composite view in
FIJI / napari. When no `markers` are provided the first 10 channels (by
index) are displayed.

## Usage

``` r
plot_qptiff(
  x,
  markers = NULL,
  resolution = 1L,
  palette = "magma",
  assay.type = "bgnorm",
  scale = c("marker", "sample", "none"),
  point_size = 1,
  pixels = c(1024L, 1024L),
  flip_y = TRUE,
  large_data_threshold = 10000L,
  ncol = NULL
)
```

## Arguments

- x:

  A [`QPTIFFImage`](QPTIFFImage.md) *or* a
  [`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html).

- markers:

  Character vector of channel names to display, or `NULL` to use the
  first 10 channels (for `QPTIFFImage`) or all features (for
  `SpatialExperiment`). At most 10 channels are shown for `QPTIFFImage`
  input; excess channels are dropped with a warning.

- resolution:

  Positive integer; pixel downsample factor for `QPTIFFImage` input
  (ignored for `SpatialExperiment`).

- palette:

  Viridis palette name used for `SpatialExperiment` input only. Default
  `"magma"`.

- assay.type:

  Assay to visualise when `x` is a `SpatialExperiment`. Default
  `"bgnorm"`.

- scale:

  Character; intensity scaling for `QPTIFFImage` input (ignored for
  `SpatialExperiment`). One of:

  `"marker"`

  :   (default) Per-channel min and 99.9th-percentile; each marker is
      stretched to full brightness independently.

  `"sample"`

  :   Global min and 99.9th-percentile computed across all channels in
      the image; preserves relative intensities between markers.

  `"none"`

  :   No scaling; values are clamped to `[0, 1]`.

- point_size:

  Base point size for `SpatialExperiment` scatter.

- pixels:

  Integer `c(width, height)` for scattermore rasterisation. Default
  `c(1024L, 1024L)`.

- flip_y:

  Logical; reverse y-axis for `SpatialExperiment`? Default `TRUE`.

- large_data_threshold:

  Cell count above which scattermore is used automatically. Default
  `10000L`.

- ncol:

  Ignored for `QPTIFFImage` (single composite panel). Number of facet
  columns for `SpatialExperiment`.

## Value

A `ggplot` object.

## Details

For a
[`SpatialExperiment`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html):
plots each cell at its spatial coordinates, coloured by its (per-channel
scaled) intensity from the requested assay using a viridis colour scale
and dark theme. When the cell count exceeds `large_data_threshold` and
scattermore is installed, the scatter layer is rasterised automatically.

## See also

[`plot_pixel_classes`](plot_pixel_classes.md),
[`plot_distributions`](plot_distributions.md),
[`bgnorm_pixels`](bgnorm_pixels.md), [`bgnorm_sce`](bgnorm_sce.md)

## Examples

``` r
path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
img  <- read_qptiff(path)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
plot_qptiff(img, markers = c("PanCK", "CD20"))


# \donttest{
# After normalisation: plot background-adjusted intensities
res <- bgnorm_pixels(img, sample_prop = 0.1)
plot_qptiff(res, markers = c("PanCK", "CD20"), scale = "sample")

# }
```
