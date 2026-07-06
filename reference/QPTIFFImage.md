# QPTIFFImage: an in-memory or on-disk multi-channel image

An S3 class representing a `[H x W x C]` multiplex image, either eagerly
loaded as a plain 3-D array or lazily backed by a `DelayedArray` (see
[`QPTIFFArraySeed-class`](QPTIFFArraySeed-class.md)). Created by
[`read_qptiff`](read_qptiff.md) or
[`as.QPTIFFImage`](as.QPTIFFImage.md), and consumed by
[`bgnorm_pixels`](bgnorm_pixels.md), [`plot_qptiff`](plot_qptiff.md),
and [`write_qptiff`](write_qptiff.md).

## Value

Not applicable; documents the `QPTIFFImage` class structure.
