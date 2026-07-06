# Changelog

## bgnormR 0.99.0

### New Features

- Added [`bgnorm_pixels()`](../reference/bgnorm_pixels.md) for
  pixel-level background normalisation of `QPTIFFImage` objects,
  supporting parallel channel processing via `BiocParallel`.
- Added [`bgnorm_cells()`](../reference/bgnorm_cells.md) and
  [`bgnorm_sce()`](../reference/bgnorm_sce.md) for cell-level
  normalisation of `SummarizedExperiment` (including
  `SingleCellExperiment` and `SpatialExperiment`) objects.
- Added `QPTIFFImage` S3 class (3-D array `[H × W × C]`) and
  [`as.QPTIFFImage()`](../reference/as.QPTIFFImage.md) coercion generic.
- Added [`read_qptiff()`](../reference/read_qptiff.md) for reading Akoya
  PhenoCycler-Fusion QPTIFF images without Java.
- Added `QPTIFFArraySeed` / `DelayedArray` backend for out-of-core lazy
  loading of QPTIFF images.
- Added unified plotting functions:
  [`plot_qptiff()`](../reference/plot_qptiff.md),
  [`plot_pixel_classes()`](../reference/plot_pixel_classes.md),
  [`plot_distributions()`](../reference/plot_distributions.md),
  [`plot_jsd_heatmap()`](../reference/plot_jsd_heatmap.md) supporting
  both `QPTIFFImage` (raster) and `SpatialExperiment` (scatter) input.
- Added [`jsd_qc()`](../reference/jsd_qc.md) and
  [`qc_summary()`](../reference/qc_summary.md) for Jensen-Shannon
  Divergence quality control.
- Added [`log_transform()`](../reference/log_transform.md),
  [`inv_log_transform()`](../reference/inv_log_transform.md), and
  [`median_filter_3x3()`](../reference/median_filter_3x3.md)
  pre-processing utilities.
