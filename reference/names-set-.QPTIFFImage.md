# Set the channel names of a QPTIFFImage

Replaces the third-dimension names (channel names) of a `QPTIFFImage`.
Works for both eager (in-memory) and lazy (on-disk) objects.

## Usage

``` r
# S3 method for class 'QPTIFFImage'
names(x) <- value
```

## Arguments

- x:

  A `QPTIFFImage`.

- value:

  Character vector of length equal to the number of channels.

## Value

`x` with updated channel names.

## Examples

``` r
arr <- array(1:8, dim = c(2, 2, 2),
             dimnames = list(NULL, NULL, c("ch1", "ch2")))
img <- as.QPTIFFImage(arr)
names(img) <- c("DAPI", "CD3")
names(img)
#> [1] "DAPI" "CD3" 
```
