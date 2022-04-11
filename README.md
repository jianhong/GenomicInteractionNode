# GenomicInteractionNodes

The `GenomicInteractionNodes` package define genomic interaction hot regions
as interaction nodes with following attributes:

1. it must contain multiple interaction loops,

2. it regulates one or more target genes.

# Installation

You can install the package via `devtools::install_github` from `github`.
```{r}
library(devtools)
install_github("jianhong/GenomicInteractionNodes")
```

You can also try to install it via `BiocManager::install` when it is ready in Bioconductor.
```{r}
library(BiocManager)
install("GenomicInteractionNodes")
```

# Documentation

To view documentation of GenomicInteractionNodes, start R and enter:

```{r}
browseVignettes("GenomicInteractionNodes")
```