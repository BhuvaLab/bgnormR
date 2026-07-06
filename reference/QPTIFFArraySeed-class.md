# Seed class for lazy on-disk QPTIFF access

An S4 class that satisfies the DelayedArray seed contract, enabling
lazy, block-based processing of large Akoya PhenoCycler-Fusion QPTIFF
images. Each channel (TIFF page) is read from disk only when data for
that channel is actually requested by the DelayedArray framework.

Called internally by the DelayedArray framework. Reads only the
requested channels and the spatial region specified by `index` from
disk, keeping memory use proportional to the data actually needed.

## Usage

``` r
# S4 method for class 'QPTIFFArraySeed'
dim(x)

# S4 method for class 'QPTIFFArraySeed'
dimnames(x)

# S4 method for class 'QPTIFFArraySeed'
type(x)

# S4 method for class 'QPTIFFArraySeed'
extract_array(x, index)
```

## Arguments

- x:

  A `QPTIFFArraySeed`.

- index:

  List of length 3. Each element is either `NULL` (all indices for that
  dimension) or a sorted integer vector of 1-based indices.

## Value

An ordinary array with dimensions
`c(length(i_rows), length(i_cols), length(i_channels))`.

## Details

Internally the seed reads the full TIFF/BigTIFF directory chain at
construction time (fast: only header bytes, no pixel data) and stores
tile/strip layout metadata for each channel page. Pixel data is fetched
per-tile (when the image is tiled) or per-strip with DEFLATE/zlib
decompression handled by R's built-in `memDecompress`. Unsupported
compression types fall back to
[`tiff::readTIFF`](https://rdrr.io/pkg/tiff/man/readTIFF.html).

Do not construct this class directly; use
[`read_qptiff`](read_qptiff.md)`(path, lazy = TRUE)`.

## Slots

- `filepath`:

  Absolute path to the QPTIFF file.

- `.dim`:

  Integer vector `[H, W, C]`.

- `.dimnames`:

  List of dimnames; element 3 holds channel names.

- `page_layouts`:

  List of per-channel page layout objects (one per channel) each
  containing the TIFF structural metadata needed to fetch that page.

- `dtype`:

  Character, `"integer"` (raw 16-bit values 0-65535) or `"double"`
  (normalised to \[0, 1\]).

- `metadata`:

  List; rich QPI metadata as returned by `.parse_qpi_xml`.

- `level`:

  Integer; pyramid resolution level (1 = full resolution).
