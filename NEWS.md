# bgnormR 0.99.1

## Improvements
* `plot_qptiff()` and `plot_pixel_classes()` now render plot titles, subtitles, and captions in white so they are legible on the dark image background, across both the `QPTIFFImage` raster and `SpatialExperiment` scatter render paths.

# bgnormR 0.99.0

## New Features
* Added `bgnorm_pixels()` for pixel-level background normalisation of `QPTIFFImage` objects, supporting parallel channel processing via `BiocParallel`.
* Added `bgnorm_cells()` and `bgnorm_sce()` for cell-level normalisation of `SummarizedExperiment` (including `SingleCellExperiment` and `SpatialExperiment`) objects.
* Added `QPTIFFImage` S3 class (3-D array `[H × W × C]`) and `as.QPTIFFImage()` coercion generic.
* Added `read_qptiff()` for reading Akoya PhenoCycler-Fusion QPTIFF images without Java.
* Added `QPTIFFArraySeed` / `DelayedArray` backend for out-of-core lazy loading of QPTIFF images.
* Added unified plotting functions: `plot_qptiff()`, `plot_pixel_classes()`, `plot_distributions()`, `plot_jsd_heatmap()` supporting both `QPTIFFImage` (raster) and `SpatialExperiment` (scatter) input.
* Added `jsd_qc()` and `qc_summary()` for Jensen-Shannon Divergence quality control.
* Added `log_transform()`, `inv_log_transform()`, and `median_filter_3x3()` pre-processing utilities.
