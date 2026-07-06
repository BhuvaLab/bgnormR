# Replace the dimnames of a QPTIFFImage

Low-level setter used by `names<-.QPTIFFImage` and direct
`dimnames(img) <- ` assignments. Preserves the `QPTIFFImage` class and
the `"metadata"` attribute.

## Usage

``` r
# S3 method for class 'QPTIFFImage'
dimnames(x) <- value
```

## Arguments

- x:

  A `QPTIFFImage`.

- value:

  A list of length 3 (`NULL` elements are allowed for the row and column
  dimensions).

## Value

`x` with updated dimnames.
