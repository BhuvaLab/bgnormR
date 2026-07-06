# Per-channel bgnorm results stored in a QPTIFFImage

Returns the named list of [`BgnormResult`](BgnormResult.md) objects
attached to a [`QPTIFFImage`](QPTIFFImage.md) by
[`bgnorm_pixels`](bgnorm_pixels.md), or `NULL` if the image has not been
background-normalised.

## Usage

``` r
bgnorm_results(x, ...)
```

## Arguments

- x:

  A `QPTIFFImage`.

- ...:

  Unused.

## Value

A named list of `BgnormResult` objects (one per processed channel), or
`NULL`.

## See also

[`bgnorm_pixels`](bgnorm_pixels.md)

## Examples

``` r
path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
img <- read_qptiff(path)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
adj <- bgnorm_pixels(img, markers = names(img)[1])
bgnorm_results(adj)
#> $CD20
#> BgnormResult (pixel-level)
#>   n = 440000 
#>   Component means: 0.027 0.385 1.151 
#>   JSD (QC metric): 0.4764 
#>   BIC (G=2): -0.35   BIC (G=3): -0.02 
#>   No signal detected: FALSE 
#>   Quantile normalised: FALSE 
#>   Tissue positivity: 18.5% 
#> 
#> $CD3e
#> BgnormResult (pixel-level)
#>   n = 440000 
#>   Component means: 0.145 0.869 1.965 
#>   JSD (QC metric): 0.4459 
#>   BIC (G=2): -1.63   BIC (G=3): -1.41 
#>   No signal detected: FALSE 
#>   Quantile normalised: FALSE 
#>   Tissue positivity: 16.3% 
#> 
#> $CD8
#> BgnormResult (pixel-level)
#>   n = 440000 
#>   Component means: 0.023 0.122 0.803 
#>   JSD (QC metric): 0.7098 
#>   BIC (G=2): 1.79   BIC (G=3): 2.39 
#>   No signal detected: FALSE 
#>   Quantile normalised: FALSE 
#>   Tissue positivity: 11.6% 
#> 
#> $PanCK
#> BgnormResult (pixel-level)
#>   n = 440000 
#>   Component means: 0.03 0.162 1.878 
#>   JSD (QC metric): 0.8813 
#>   BIC (G=2): -1.69   BIC (G=3): -1.4 
#>   No signal detected: FALSE 
#>   Quantile normalised: FALSE 
#>   Tissue positivity: 67.9% 
#> 
#> $Vimentin
#> BgnormResult (pixel-level)
#>   n = 440000 
#>   Component means: 0.015 0.318 2.013 
#>   JSD (QC metric): 0.7546 
#>   BIC (G=2): -2.22   BIC (G=3): -1.88 
#>   No signal detected: FALSE 
#>   Quantile normalised: FALSE 
#>   Tissue positivity: 49.6% 
#> 
```
