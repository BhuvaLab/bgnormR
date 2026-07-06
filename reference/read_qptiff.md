# Read an Akoya PhenoCycler-Fusion QPTIFF image

Reads a QPTIFF file produced by the Akoya PhenoCycler-Fusion (formerly
CODEX) platform, Cell DIVE, or Vectra / Polaris scanners into R without
Java. Handles three live format variants: brightfield RGB, Polaris
ScanBand XML, and Fusion per-page JSON+XML.

## Usage

``` r
read_qptiff(path, channels = NULL, level = 1L, as_integer = TRUE, lazy = FALSE)
```

## Arguments

- path:

  Character, path to the QPTIFF file.

- channels:

  Character vector of channel names to load, an integer vector of
  1-based channel indices, or `NULL` (default) to load all channels.

- level:

  Integer, pyramid resolution level. `1` = full resolution (default),
  `2` = half resolution, etc.

- as_integer:

  Logical; return raw 16-bit integers (0-65535) rather than normalised
  `[0, 1]` doubles? Default `TRUE`.

- lazy:

  Logical; if `TRUE` return a
  [`DelayedArray`](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html)
  backed by a [`QPTIFFArraySeed`](QPTIFFArraySeed-class.md) that reads
  pages from disk on demand. If `FALSE` (default) load all requested
  channels into memory and return a [`QPTIFFImage`](QPTIFFImage.md).

## Value

- `lazy = FALSE`:

  A `QPTIFFImage` - a 3-D numeric array `[height, width, channels]` with
  class `c("QPTIFFImage", "array")`. Channel names are stored in
  `dimnames(img)[[3]]`. Standard array subscripting works:
  `img[, , "DAPI"]` extracts a single channel as a matrix;
  `img[1:512, 1:512, ]` crops spatially. Rich metadata is in
  `attr(img, "metadata")`.

- `lazy = TRUE`:

  A `DelayedArray` with `dim = c(H, W, C)`. Individual channel pages are
  read from disk only when accessed. The seed is a `QPTIFFArraySeed`
  accessible via `DelayedArray::seed(arr)`.

## Details

Channel names and rich per-channel metadata (fluorophore, exposure time,
wavelengths, filters) are parsed from the TIFF ImageDescription tag
(270) of each page. The PerkinElmer XML schema is detected
automatically.

## Examples

``` r
path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
img  <- read_qptiff(path)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
dim(img)                      # c(H, W, n_channels)
#> [1] 550 800   5
names(img)                    # channel names
#> [1] "CD20"     "CD3e"     "CD8"      "PanCK"    "Vimentin"
cd20 <- img[, , "CD20"]      # extract one channel as a 2-D matrix

# Load specific channels only
img2 <- read_qptiff(path, channels = c("CD20", "PanCK"))
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 2 channel(s) ...
names(img2)
#> [1] "CD20"  "PanCK"

# \donttest{
# Lazy / out-of-core load (backed by DelayedArray)
arr <- read_qptiff(path, lazy = TRUE)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
# }
```
