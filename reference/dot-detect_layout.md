# Determine number of channels and pyramid levels from IFD descriptions

Detects the repeating-page pattern:

- All identical -\> Polaris; channel count from ScanBands-i.

- Repeating with period N -\> Fusion paged; N channels.

## Usage

``` r
.detect_layout(descriptions, samples_per_pixel)
```
