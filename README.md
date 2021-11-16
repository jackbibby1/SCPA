
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCPA

<!-- badges: start -->
<!-- badges: end -->

![SCPA outline](images/scpa_outline.png)

SCPA is a method for pathway analysis in single cell RNA-seq data. It’s
based on a novel approach to pathway analysis that defines pathway
activity as a change in multivariate distribution of a given pathway
across conditions, rather than enrichment or overrepresentation of
genes. This approach allows for a number of benefits over current
methods, but mainly:

1.  You will identify pathways that show enrichment in a given
    population AND also identify pathways with no overall enrichment but
    alterations in the multivariate distribution of that pathway. You
    essentially get the best of both worlds, as pathways with changes in
    multivariate distribution but no overall enrichment are still
    interestingly different pathways.

2.  You can compare multiple conditions simultaneously e.g. compare
    across 3 time points, or across multiple phases of a pseuodotime
    trajectory. This means you can assess pathway activity through
    multiple stimulation times, or across cell differentiation

## 1. Installation

You can install the development version of SCPA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jackbibby1/SCPA")
```

## 2. Tutorials
