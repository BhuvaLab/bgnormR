# bgnormR

**bgnormR** implements *bgnorm* (Kharbanda, Tubelleza *et al.* 2025), a
generative statistical framework for background correction,
normalisation, and quality control of multiplex spatial proteomics data
from Akoya PhenoCycler-Fusion, Cell DIVE, IMC, and CosMx platforms. The
package includes a **Java-free** reader and writer for Akoya QPTIFF
images and provides full integration with Bioconductor’s
[SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/)
/
[SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment/)
ecosystem.

## How it works

bgnorm models log₂-transformed fluorescence intensities as a
**three-component Gaussian Mixture Model** (GMM) representing
background, non-specific binding, and biological signal. A closed-form
deconvolution step isolates the signal component, and the Jensen-Shannon
Divergence (JSD) between components 2 and 3 provides an automated
staining-quality metric.

![Three-component signal model](articles/images/bgnorm.png)

Three-component signal model

## Installation

Install from Bioconductor (once accepted):

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("bgnormR")
```

Install the development version from GitHub:

``` r

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("BhuvaLab/bgnormR")
```

## Quick start

### Pixel-level normalisation

``` r

library(bgnormR)

# Read a multi-channel QPTIFF (no Java required)
img <- read_qptiff("path/to/image.qptiff")
dim(img)    # height × width × channels
names(img)  # protein panel

# Fit 3-component GMM and background-correct every channel
res <- bgnorm_pixels(img, sample_prop = 0.1)

# Inspect the model for one channel
bgnorm_results(res)[["PanCK"]]

# Visualise distributions and class assignments
plot_distributions(res)
plot_pixel_classes(res, markers = c("PanCK", "CD20"))

# Quality control: JSD per marker
qc_summary(res)
plot_jsd_heatmap(res)

# Export the corrected image
write_qptiff(res, "path/to/output.qptiff")
```

### Cell-level normalisation (SummarizedExperiment)

``` r

library(SummarizedExperiment)

# se: a SummarizedExperiment with raw counts assay
se <- bgnorm_sce(se, assay.type = "counts", name = "bgnorm")

# Per-marker model parameters are stored in metadata
metadata(se)$bgnorm_results[["CD20"]]$jsd
```

## Key functions

| Function | Description |
|----|----|
| [`read_qptiff()`](reference/read_qptiff.md) | Read Akoya QPTIFF images (eager or lazy/DelayedArray) |
| [`write_qptiff()`](reference/write_qptiff.md) | Write a `QPTIFFImage` to multi-page 16-bit TIFF with embedded metadata |
| [`bgnorm_pixels()`](reference/bgnorm_pixels.md) | Pixel-level 3-component GMM correction on a `QPTIFFImage` |
| [`bgnorm_cells()`](reference/bgnorm_cells.md) | Cell-level 2-component GMM correction on a numeric vector |
| [`bgnorm_sce()`](reference/bgnorm_sce.md) | Apply `bgnorm_cells` to every marker in a `SummarizedExperiment` |
| [`qc_summary()`](reference/qc_summary.md) | Per-marker JSD and signal proportion table |
| [`plot_distributions()`](reference/plot_distributions.md) | Histogram + fitted GMM density curves per marker |
| [`plot_pixel_classes()`](reference/plot_pixel_classes.md) | Spatial class-assignment map (Background / Non-specific / Signal) |
| [`plot_qptiff()`](reference/plot_qptiff.md) | Multi-channel intensity composite image |
| [`plot_jsd_heatmap()`](reference/plot_jsd_heatmap.md) | JSD quality-control heatmap across markers and samples |

## Citation

Kharbanda M, Tubelleza R, Salim A, Bhuva DD (2025). *bgnorm: a
generative framework for background correction and normalisation of
multiplex spatial proteomics.* bioRxiv.
<https://doi.org/10.1101/2025.placeholder>

## License

MIT © Dharmesh D. Bhuva
