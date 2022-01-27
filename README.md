
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

# Single Cell Pathway Analysis

![](man/figures/scpa_outline.png)

### About SCPA

SCPA is a method for pathway analysis in single cell RNA-seq data. It’s
based on a novel approach to pathway analysis that defines pathway
activity as a change in multivariate distribution of a given pathway
across conditions, rather than enrichment or over representation of
genes. This approach allows for a number of benefits over current
methods, including:

1.  You will identify pathways that show enrichment in a given
    population AND also identify pathways with no overall enrichment but
    alterations in the multivariate distribution of that pathway. You
    essentially get the best of both worlds, as pathways with changes in
    multivariate distribution but no overall enrichment are still
    interestingly different pathways.

2.  SCPA allows for multisample testing, so you can compare multiple
    conditions simultaneously e.g. compare across 3 time points, or
    across multiple phases of a [pseuodotime trajectory](). This means
    you can assess pathway activity through multiple stimulation phases,
    or across cell differentiation

To see the stats behind SCPA, you can see our paper in JASA
[here](https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1791131)

We are currently putting a paper together to show the power of SCPA
applied in a systematic way over a T cell scRNAseq dataset, which will
be available on bioRxiv soon.

## 1. Installation

You can install SCPA by running:

``` r
# install.packages("devtools")
devtools::install_github("jackbibby1/SCPA")
```

## 2. Tutorials

For tutorials on how to use SCPA, [click this
link](jackbibby1.github.io/SCPA)

## 3. Issues

To report any issues or submit comments please use:
<https://github.com/jackbibby1/SCPA/issues>
