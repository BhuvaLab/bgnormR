# Parse a PerkinElmer QPI ImageDescription XML

Parse a PerkinElmer QPI ImageDescription XML

## Usage

``` r
.parse_qpi_xml(
  xml_str,
  per_page_xmls = NULL,
  n_channels = NULL,
  is_rgb = FALSE
)
```

## Arguments

- xml_str:

  Raw XML string from TIFF tag 270 of the first page.

- per_page_xmls:

  Character vector of XML strings, one per channel page.

- n_channels:

  Fallback channel count when XML is absent.

- is_rgb:

  TRUE when the first page is RGB (brightfield signal).

## Value

Named list with `slide`, `image_info`, `channels`, `format`, `raw_xml`.
