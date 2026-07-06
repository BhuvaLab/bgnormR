# Fallback page reader via tiff::readTIFF

Used when compression is not natively supported (LZW, JPEG, etc.). Loads
ALL pages into memory then extracts the requested page; this is
memory-intensive for large multi-channel files.

## Usage

``` r
.read_page_fallback(path, layout, rows, cols)
```
