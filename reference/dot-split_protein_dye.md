# Split a compound biomarker name into protein and fluorescence dye

Parses names of the form `"CD3e-AF647"` into a protein component
(`"CD3e"`) and a fluorescence dye component (`"AF647"`). The
already-parsed `fluorophore` field is used as the primary signal; a
regex against known dye naming conventions serves as a fallback.

## Usage

``` r
.split_protein_dye(name, fluorophore = NULL)
```

## Arguments

- name:

  Character scalar; raw channel / biomarker name.

- fluorophore:

  Character scalar or `NULL`; fluorophore already parsed from the XML
  `<Fluorophore>` element.

## Value

A list with elements `protein` (character) and `dye` (character or
`NULL`).

## Details

Handles multi-hyphen names correctly: `"HLA-A-Atto550"` -\> protein
`"HLA-A"`, dye `"Atto550"`.
