# Access the metadata embedded in a QPTIFFImage

Retrieves the named list of slide / image / channel metadata stored
inside a [`QPTIFFImage`](QPTIFFImage.md). The generic falls back to
[`metadata`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html) for
S4 objects such as `SummarizedExperiment`, so `metadata(img)` and
`metadata(spe)` both work when bgnormR is loaded.

## Usage

``` r
metadata(x, ...)

# S3 method for class 'QPTIFFImage'
metadata(x, ...)
```

## Arguments

- x:

  A [`QPTIFFImage`](QPTIFFImage.md).

- ...:

  Unused.

## Value

A named list of metadata.

## Details

The returned list typically contains:

- `format`:

  One of `"brightfield"`, `"polaris_scanband"`, or `"fusion_paged"`.

- `channels`:

  A list of per-channel lists, each with fields `name` (protein name,
  dye suffix removed), `fluorophore` (parsed from the `<Fluorophore>`
  XML tag), `dye_from_name` (the dye suffix that was stripped from the
  compound biomarker name), `exposure_time_us`,
  `emission_wavelength_nm`, etc.

- `image_info`:

  Physical image properties (pixel size, etc.).

- `n_levels`:

  Number of pyramid resolution levels.

## Examples

``` r
path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
img  <- read_qptiff(path)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
meta <- metadata(img)
#> Error: unable to find an inherited method for function ‘metadata’ for signature ‘x = "QPTIFFImage"’
meta$format
#> Error: object 'meta' not found
# Per-channel dye information
ch_meta <- meta$channels
#> Error: object 'meta' not found
data.frame(
  protein = sapply(ch_meta, `[[`, "name"),
  dye     = sapply(ch_meta, function(ch) ch$dye_from_name %||% NA_character_)
)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'X' in selecting a method for function 'sapply': object 'ch_meta' not found
```
