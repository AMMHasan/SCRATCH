# scRatch: Copy number change for Relationship Assessment and Testing Clone Histories

<!-- badges: start -->
<!-- badges: end -->

scRatch is a an implementation of hierarchical clustering based algorithm to determine the evolutionary relationship of metastases based on correlation of the copy number at the transition points, which are the genomic positions of copy number changes between adjacent genomic segments derived from copy number analysis on whole genome sequencing data. This is a beta version before an official release.

## Installation

You can install the development version of scRatch from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AMMHasan/scRatch")
```

## Example

This is a basic example which shows how to create a metastatic relationship plot from the bam files derived from whole genome sequencing of tumours:

``` r
library(scRatch)
## basic example code

BAM_path <- "~/Documents/Research/prostate_cancer/PEACE/samples/tumour_BAM_subset/"
bin_size <- 500
pattern <- "PEA310"
seed <- 101


CN_called_obj <- 
  generate_copyNumberCalled_obj(BAM_path, pattern, bin_size, seed)

TP_CN_mat <- 
  generate_TP_CN_matrix(BAM_path, pattern, CN_called_obj) 

build_tree(TP_CN_mat) %>% 
  plot
```



