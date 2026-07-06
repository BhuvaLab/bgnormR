# Lazy constructor for QPTIFFImage (DelayedArray backing)

Lazy constructor for QPTIFFImage (DelayedArray backing)

## Usage

``` r
.new_QPTIFFImage_lazy(da, metadata = list(), bgnorm_results = NULL)
```

## Arguments

- da:

  A `DelayedArray` from a `QPTIFFArraySeed`.

- metadata:

  List of slide / image / channel metadata.

- bgnorm_results:

  Named list of `BgnormResult` objects, or `NULL`.
