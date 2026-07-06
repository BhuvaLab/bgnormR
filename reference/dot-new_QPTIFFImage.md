# Eager constructor for QPTIFFImage (in-memory array backing)

Eager constructor for QPTIFFImage (in-memory array backing)

## Usage

``` r
.new_QPTIFFImage(arr, metadata = list(), bgnorm_results = NULL)
```

## Arguments

- arr:

  Numeric array of dimension `c(H, W, C)`.

- metadata:

  List of slide / image / channel metadata.

- bgnorm_results:

  Named list of `BgnormResult` objects (one per channel), or `NULL`.
