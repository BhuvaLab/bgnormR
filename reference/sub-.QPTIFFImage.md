# Subset a QPTIFFImage

Standard 3-D array subscripting with all missing-subscript combinations
supported. The result is *always* a `QPTIFFImage`, including when a
single channel is selected. To extract a plain 2-D matrix for one
channel, use `unclass(img[, , "DAPI"])[, , 1L]`.

- `img[, , "DAPI"]` - single channel (returns a 1-channel QPTIFFImage)

- `img[, , c("DAPI","CD3")]` - multi-channel sub-image

- `img[1:512, 1:512, ]` - spatial crop

## Usage

``` r
# S3 method for class 'QPTIFFImage'
x[i, j, k, drop = FALSE]
```

## Arguments

- x:

  A `QPTIFFImage`.

- i:

  Row indices (spatial), or missing.

- j:

  Column indices (spatial), or missing.

- k:

  Channel indices (numeric or character), or missing.

- drop:

  Ignored; the result is always a 3-D `QPTIFFImage`.

## Value

A `QPTIFFImage`.
