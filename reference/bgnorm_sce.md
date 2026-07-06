# Apply cell-level bgnorm to a SingleCellExperiment or matrix

Runs [`bgnorm_cells`](bgnorm_cells.md) independently for each marker
(column) and returns an updated object with adjusted intensities in a
new assay.

## Usage

``` r
bgnorm_sce(
  x,
  assay.type = "counts",
  name = "bgnorm",
  cofactor = 150,
  quantile_norm = FALSE,
  quantile = 0.75,
  BPPARAM = BiocParallel::SerialParam(),
  ...
)
```

## Arguments

- x:

  A `SummarizedExperiment` (including `SingleCellExperiment` and
  `SpatialExperiment` subclasses) or a numeric matrix (cells \\\times\\
  markers).

- assay.type:

  Character; name of the assay to normalise. Default `"counts"`.

- name:

  Character; name for the output assay (SCE) or returned matrix
  attribute. Default `"bgnorm"`.

- cofactor:

  Cofactor for log2 transform (default 150).

- quantile_norm:

  Logical; apply bgnormQ?

- quantile:

  Quantile for normalisation.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance.

- ...:

  Additional arguments forwarded to `bgnorm_cells`.

## Value

The input object with a new assay (`name`) holding the adjusted
intensities. The per-marker `BgnormResult` list is stored in
`metadata(x)$bgnorm_results`.

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:bgnormR’:
#> 
#>     metadata
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
set.seed(5)
counts <- matrix(
  exp(rnorm(2000, log(300), 0.9)),
  nrow = 10, ncol = 200,
  dimnames = list(paste0("marker", seq_len(10)), paste0("cell", seq_len(200)))
)
se <- SummarizedExperiment(assays = list(counts = counts))
se <- bgnorm_sce(se, assay.type = "counts", name = "bgnorm")
assayNames(se)
#> [1] "counts" "bgnorm"
```
