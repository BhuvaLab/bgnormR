# BgnormResult S3 class

A lightweight list-based S3 class returned by
[`bgnorm_pixels`](bgnorm_pixels.md) (per channel, via
[`bgnorm_results`](bgnorm_results.md)) and
[`bgnorm_cells`](bgnorm_cells.md).

## Value

Not applicable; documents the `BgnormResult` class structure.

## Fields

- `parameters`:

  List with elements `means`, `sds`, `props`.

- `n`:

  Integer; total number of input observations (pixels or cells,
  including zeros) for this channel.

- `threshold`:

  For the 3-component (pixel-level) model: a single numeric scalar
  giving the maximum adjusted intensity of Non-specific pixels. Combined
  with the implicit class-1 boundary at 0, this encodes the full
  classification: `adj == 0` -\> Background; `0 < adj <= threshold` -\>
  Non-specific; `adj > threshold` -\> Signal. `NULL` for 2-component
  models (cell-level or no-signal fallback), where the only boundary is
  0 (adjusted \> 0 -\> Signal).

- `jsd`:

  Jensen-Shannon Divergence QC metric between components 2 and 3.
  `NA_real_` for no-signal channels.

- `level`:

  `"pixel"` or `"cell"`.

- `quantile_norm`:

  Logical; whether quantile normalisation was applied.

- `histogram`:

  List with `$breaks` and `$density` vectors from a pre-computed
  histogram of log2-transformed intensities. Used by
  [`plot_distributions`](plot_distributions.md).

- `no_signal`:

  Logical; `TRUE` when the 3-component GMM failed and the channel was
  fitted with a 2-component fallback, indicating no detectable
  biological signal.

- `bic`:

  Named numeric vector `c(G2 = ..., G3 = ...)` with mclust BIC values
  (higher is better in mclust's convention) divided by the number of
  points used to fit, giving a point-specific estimate comparable across
  channels regardless of sample size. `NULL` for cell-level results.
