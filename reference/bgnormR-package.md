# bgnormR: Background Normalisation for Multiplex Spatial Proteomics

Implements *bgnorm*, a generative statistical framework for background
correction, normalisation, and quality control in multiplex spatial
proteomics (Kharbanda et al. 2025).

## Details

The framework models log-transformed fluorescence intensities as a
three-component Gaussian mixture representing:

1.  Background (empty space / instrument noise)

2.  Non-specific binding and autofluorescence

3.  True biological signal

Core functions:

- [`bgnorm_pixels`](bgnorm_pixels.md) - pixel-level background
  correction

- [`bgnorm_cells`](bgnorm_cells.md) - cell-level background correction

- [`bgnorm_sce`](bgnorm_sce.md) - normalise a SummarizedExperiment /
  SpatialExperiment

- [`jsd_qc`](jsd_qc.md) - Jensen-Shannon Divergence quality metric

- [`read_qptiff`](read_qptiff.md) - Akoya PhenoCycler-Fusion image
  reader

- [`plot_qptiff`](plot_qptiff.md),
  [`plot_pixel_classes`](plot_pixel_classes.md),
  [`plot_distributions`](plot_distributions.md),
  [`plot_jsd_heatmap`](plot_jsd_heatmap.md) - plotting utilities

## See also

Useful links:

- <https://github.com/BhuvaLab/bgnormR>

- Report bugs at <https://github.com/BhuvaLab/bgnormR/issues>

## Author

**Maintainer**: Dharmesh D. Bhuva <dhrmsh36@gmail.com>
([ORCID](https://orcid.org/0000-0002-6398-9157))

Authors:

- Dharmesh D. Bhuva <dhrmsh36@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6398-9157))

- Rafael Tubelleza

- Malvika Kharbanda

- Agus Salim
