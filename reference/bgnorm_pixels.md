# Pixel-level background normalisation (bgnorm) for a QPTIFFImage

Fits a three-component Gaussian Mixture Model (GMM) to the
log2-transformed pixel intensities for each channel in a
[`QPTIFFImage`](QPTIFFImage.md) and deconvolves the biological signal
component from background sources. All channels are processed
independently; use `channels` to restrict processing to a subset.
Optionally runs channels in parallel via BiocParallel.

## Usage

``` r
bgnorm_pixels(
  img,
  channels = NULL,
  cofactor = 150,
  quantile_norm = FALSE,
  quantile = 0.75,
  sample_prop = 0.1,
  BPPARAM = BiocParallel::SerialParam(),
  ...
)
```

## Arguments

- img:

  A [`QPTIFFImage`](QPTIFFImage.md) `[H x W x C]`. Use
  [`as.QPTIFFImage`](as.QPTIFFImage.md) to promote a plain 3-D array or
  2-D matrix.

- channels:

  Character or integer vector of channels to process, or `NULL`
  (default) to process all channels.

- cofactor:

  Positive numeric cofactor for the log2 transform. Default 150 is
  appropriate for 16-bit Akoya PhenoCycler-Fusion images.

- quantile_norm:

  Logical; apply model-based quantile normalisation (bgnormQ)? Default
  `FALSE`.

- quantile:

  Quantile of the signal component used for normalisation when
  `quantile_norm = TRUE`. Default 0.75.

- sample_prop:

  Numeric in (0, 1\]; proportion of non-zero pixels used to fit the GMM.
  Values \< 1 draw a random subset, accelerating fitting on large
  images. Default `1` (all non-zero pixels).

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance. Default
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html).

- ...:

  Additional arguments forwarded to the internal GMM fitter.

## Value

A `QPTIFFImage` `[H x W x C]` of background-adjusted log-intensities.
The per-channel [`BgnormResult`](BgnormResult.md) objects (GMM
parameters, JSD) are stored as an attribute; retrieve them with
[`bgnorm_results`](bgnorm_results.md)`(result)`.

## Details

The preprocessing, model, and deconvolution steps follow Kharbanda et
al. 2025:

1.  Log2-transform: \\I\_{\log} = \log_2(I/c + 1)\\.

2.  Fit a 3-component GMM via EM (mclust), filtering zeros and
    optionally subsampling for speed.

3.  Resolve component identities (background, non-specific, signal).

4.  Apply the deconvolution formula to recover \\U_3\\.

5.  Optionally normalise by the 75th-percentile of the signal component
    (bgnormQ).

## See also

[`bgnorm_cells`](bgnorm_cells.md),
[`bgnorm_results`](bgnorm_results.md), [`jsd_qc`](jsd_qc.md),
[`as.QPTIFFImage`](as.QPTIFFImage.md)

## Examples

``` r
set.seed(42)
px  <- c(exp(rnorm(500, log(50),   0.4)),
         exp(rnorm(400, log(300),  0.6)),
         exp(rnorm(100, log(2000), 1.0)))
arr <- array(c(px, px), dim = c(length(px), 1L, 2L),
             dimnames = list(NULL, NULL, c("DAPI", "CD3")))
img <- as.QPTIFFImage(arr)
res <- bgnorm_pixels(img)
#> Sampling not required: using all 1000 non-zero pixels for GMM fitting.
#> Sampling not required: using all 1000 non-zero pixels for GMM fitting.
names(res)                        # channel names of the adjusted image
#> [1] "DAPI" "CD3" 
bgnorm_results(res)[["DAPI"]]     # BgnormResult for the DAPI channel
#> BgnormResult (pixel-level)
#>   n = 1000 
#>   Component means: 0.415 1.516 4.248 
#>   JSD (QC metric): 0.722 
#>   BIC (G=2): -2.34   BIC (G=3): -2.04 
#>   No signal detected: FALSE 
#>   Quantile normalised: FALSE 
#>   Tissue positivity: 17.4% 
# Process only selected channels
res2 <- bgnorm_pixels(img, channels = "DAPI")
#> Sampling not required: using all 1000 non-zero pixels for GMM fitting.

# \donttest{
# Real tissue image shipped with the package
path <- system.file("extdata", "PA_HNC_sample.qptiff", package = "bgnormR")
img  <- read_qptiff(path)
#> Reading TIFF directory structure ...
#> Reading IFD page layouts ...
#> Loading 5 channel(s) ...
res  <- bgnorm_pixels(img, sample_prop = 0.1)
bgnorm_results(res)[["PanCK"]]    # per-channel model summary
#> BgnormResult (pixel-level)
#>   n = 440000 
#>   Component means: 0.03 0.161 1.886 
#>   JSD (QC metric): 0.8818 
#>   BIC (G=2): -1.67   BIC (G=3): -1.38 
#>   No signal detected: FALSE 
#>   Quantile normalised: FALSE 
#>   Tissue positivity: 67% 
# }
```
