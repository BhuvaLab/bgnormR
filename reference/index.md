# Package index

## All functions

- [`BgnormResult`](BgnormResult.md) : BgnormResult S3 class
- [`dim(`*`<QPTIFFArraySeed>`*`)`](QPTIFFArraySeed-class.md)
  [`dimnames(`*`<QPTIFFArraySeed>`*`)`](QPTIFFArraySeed-class.md)
  [`type(`*`<QPTIFFArraySeed>`*`)`](QPTIFFArraySeed-class.md)
  [`extract_array(`*`<QPTIFFArraySeed>`*`)`](QPTIFFArraySeed-class.md) :
  Seed class for lazy on-disk QPTIFF access
- [`QPTIFFImage`](QPTIFFImage.md) : QPTIFFImage: an in-memory or on-disk
  multi-channel image
- [`as.QPTIFFImage()`](as.QPTIFFImage.md) : Coerce a 3-D array to a
  QPTIFFImage
- [`as.array(`*`<QPTIFFImage>`*`)`](as.array.QPTIFFImage.md) :
  Materialise a QPTIFFImage to a plain 3-D array
- [`as.list(`*`<QPTIFFImage>`*`)`](as.list.QPTIFFImage.md) : Convert a
  QPTIFFImage to a named list of 2-D channel matrices
- [`as.matrix(`*`<QPTIFFImage>`*`)`](as.matrix.QPTIFFImage.md) : Coerce
  a single-channel QPTIFFImage to a matrix
- [`bgnormR-package`](bgnormR-package.md)
  [`bgnormR`](bgnormR-package.md) : bgnormR: Background Normalisation
  for Multiplex Spatial Proteomics
- [`bgnorm_cells()`](bgnorm_cells.md) : Cell-level background
  normalisation (bgnorm)
- [`bgnorm_pixels()`](bgnorm_pixels.md) : Pixel-level background
  normalisation (bgnorm) for a QPTIFFImage
- [`bgnorm_results()`](bgnorm_results.md) : Per-channel bgnorm results
  stored in a QPTIFFImage
- [`bgnorm_sce()`](bgnorm_sce.md) : Apply cell-level bgnorm to a
  SingleCellExperiment or matrix
- [`` `dimnames<-`( ``*`<QPTIFFImage>`*`)`](dimnames-set-.QPTIFFImage.md)
  : Replace the dimnames of a QPTIFFImage
- [`inv_log_transform()`](inv_log_transform.md) : Inverse log2 transform
  with cofactor
- [`jsd_gaussians()`](jsd_gaussians.md) : Jensen-Shannon Divergence
  between two Gaussian distributions
- [`jsd_qc()`](jsd_qc.md) : Compute the JSD quality control metric from
  bgnorm model parameters
- [`length(`*`<QPTIFFImage>`*`)`](length.QPTIFFImage.md) : Number of
  channels in a QPTIFFImage
- [`log_transform()`](log_transform.md) : Log2 transform with cofactor
- [`median_filter_3x3()`](median_filter_3x3.md) : Apply a 3 x 3 median
  filter to a matrix
- [`metadata()`](metadata.md) : Access the metadata embedded in a
  QPTIFFImage
- [`` `names<-`( ``*`<QPTIFFImage>`*`)`](names-set-.QPTIFFImage.md) :
  Set the channel names of a QPTIFFImage
- [`names(`*`<QPTIFFImage>`*`)`](names.QPTIFFImage.md) : Channel names
  of a QPTIFFImage (third-dimension names of the array)
- [`plot_distributions()`](plot_distributions.md) : Per-marker intensity
  distribution with fitted GMM densities
- [`plot_jsd_heatmap()`](plot_jsd_heatmap.md) : QC heatmap of
  Jensen-Shannon Divergence across markers and samples
- [`plot_pixel_classes()`](plot_pixel_classes.md) : Spatial GMM class
  assignment plot for a QPTIFFImage or SpatialExperiment
- [`plot_qptiff()`](plot_qptiff.md) : Spatial intensity plot for a
  QPTIFFImage or SpatialExperiment
- [`qc_summary()`](qc_summary.md) : Summarise QC across multiple bgnorm
  results
- [`read_qptiff()`](read_qptiff.md) : Read an Akoya PhenoCycler-Fusion
  QPTIFF image
- [`` `[`( ``*`<QPTIFFImage>`*`)`](sub-.QPTIFFImage.md) : Subset a
  QPTIFFImage
- [`write_qptiff()`](write_qptiff.md) : Write a QPTIFFImage to a TIFF
  file
