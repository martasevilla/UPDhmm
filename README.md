
# UPDhmm

<!-- badges: start -->
<!-- badges: end -->

We have developed `UPDhmm` R/Bioconductor package. The package provides a tool method to detect, classify and stablish the location of uniparental disomy events.

## Installation

You can install the current release version of UPDhmm from [Bioconductor](https://bioconductor.org/) with:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("UPDhmm")
```

You can install the development version of UPDhmm from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("martasevilla/UPDhmm")
devtools::install_github("martasevilla/UPDhmm")
```

## Example

This is a example which show you how to solve use the package in a trio-vcf file:

``` r
library(UPDhmm)
library(VariantAnnotation)
file <- system.file(package = "UPDhmm", "extdata", "test_het_mat.vcf.gz")
vcf <- readVcf(file)
processedVcf <- vcfCheck(
vcf,
proband = "NA19675",
mother = "NA19678",
father = "NA19679"
)
updEvents <- calculateEvents(largecollapsed=processedVcf)
```

