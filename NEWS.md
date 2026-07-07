# bgnormR 0.99.3

## New Features
* `read_qptiff()` now auto-detects and reads **OME-TIFF** and **OME-Zarr** (OME-NGFF) images in addition to QPTIFF. OME-TIFF channel names are read from the `<Channel Name="...">` attributes of the OME-XML; OME-Zarr stores are read (lazily, via a new `OMEZarrArraySeed` DelayedArray backend) using the optional `Rarr` package, with channel names and colours taken from `omero.channels`.
* Image metadata is now returned as a `QPTIFFMetadata` object - an OME-organised hierarchy (`slide` / `images[]` -> `image_info` / `channels[]` / `scales[]`) - with accessors `qpi_format()`, `qpi_channels()`, `qpi_channel_names()`, `qpi_pixel_size_um()`, `qpi_n_levels()`, `qpi_is_brightfield()`, and a tidy `channel_table()`.
* `write_qptiff()` now emits richer OME-TIFF metadata - physical pixel size, per-channel fluorophore / colour / emission & excitation wavelengths, per-channel exposure `<Plane>` elements, and a `qpi://vectra` `MapAnnotation` carrying PerkinElmer-specific fields - so a read -> write -> read round trip preserves this metadata. Files continue to declare `SizeC` channels (with `SizeZ = SizeT = 1`) so viewers such as QuPath read the pages as channels rather than Z/T slices.

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
