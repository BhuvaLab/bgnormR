# Write a QPTIFFImage to a TIFF file

Writes a [`QPTIFFImage`](QPTIFFImage.md) as a multi-page 16-bit
grayscale TIFF (one page per channel). Per-channel bgnorm results are
embedded as JSON in the TIFF `ImageDescription` (tag 270) of each page
so the file is self-documenting.

## Usage

``` r
write_qptiff(x, path)
```

## Arguments

- x:

  A [`QPTIFFImage`](QPTIFFImage.md) (eager or lazy).

- path:

  Output file path (character scalar). The directory must exist.

## Value

`path`, invisibly.

## Intensity transform

- bgnorm-adjusted images:

  The QPTIFFImage stores background-adjusted log\\\_2\\-intensities.
  These are inverted with \\2^x\\ before writing so the output values
  are in a linear intensity scale.

- Raw images (no bgnorm results):

  Written as-is; values are assumed to be in `[0, 65535]`.

All output values are rounded to the nearest integer and clamped to the
`[0, 65535]` range before being written as 16-bit unsigned integers.

## Metadata

Each page's `ImageDescription` tag (TIFF tag 270) is a minimal XML
document understood by [`read_qptiff`](read_qptiff.md):


    <PerkinElmerQPI>
      <Biomarker>CD20</Biomarker>
      <transform>2^x</transform>
      <bgnorm>{"level":"pixel","jsd":0.35,...}</bgnorm>
    </PerkinElmerQPI>

The `<Biomarker>` element is the channel name; `<bgnorm>` holds a JSON
object with GMM parameters, JSD, threshold, and normalisation flags
(only present for bgnorm-adjusted images). If jsonlite is not installed,
the `<bgnorm>` element is omitted.

## See also

[`read_qptiff`](read_qptiff.md), [`bgnorm_pixels`](bgnorm_pixels.md)

## Examples

``` r
path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
img  <- read_qptiff(path)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
res  <- bgnorm_pixels(img, sample_prop = 0.1)
out  <- file.path(tempdir(), "PA_HNC_bgnorm.tif")
write_qptiff(res, out)

# Round-trip: channel names and dimensions are preserved
img2 <- read_qptiff(out)
#> Reading TIFF directory structure ...
#> Note: XML root 'OME' not a recognised QPI description.
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
names(img2)
#> [1] "CD20"     "CD3e"     "CD8"      "PanCK"    "Vimentin"
dim(img2)
#> [1] 550 800   5
```
