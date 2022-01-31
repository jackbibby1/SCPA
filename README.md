
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

# Single Cell Pathway Analysis <img src="man/figures/logo.png" align="right" width=80px/>

![](man/figures/scpa_outline.png)

### About SCPA

SCPA is a method for pathway analysis in single cell RNA-seq data. It’s
a novel approach to pathway analysis that defines pathway activity as a
change in multivariate distribution of a given pathway across
conditions, rather than enrichment or over representation of genes. This
approach allows for a number of benefits over current methods,
including:

1.  You will identify pathways that show enrichment in a given
    population AND also identify pathways with no overall enrichment but
    alterations in the multivariate distribution of that pathway. You
    essentially get the best of both worlds, as pathways with changes in
    multivariate distribution but no overall enrichment are still
    interestingly different pathways.

2.  SCPA allows for multisample testing, so you can compare multiple
    conditions simultaneously e.g. compare across 3 time points, or
    across multiple phases of a [pseuodotime
    trajectory](https://jackbibby1.github.io/SCPA/articles/pseudotime.html).
    This means you can assess pathway activity through multiple
    stimulation phases, or across cell differentiation

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

If you’re viewing this page on GitHub, the SCPA webpage with all the
documentation and tutorials is
[here](https://jackbibby1.github.io/SCPA/)

We have various examples and walkthroughs, including:

-   A generic quick start
    [tutorial](https://jackbibby1.github.io/SCPA/articles/quick_start.html)
-   A
    [tutorial](https://jackbibby1.github.io/SCPA/articles/using_gene_sets.html)
    on how to get and use gene sets with SCPA
-   A
    [tutorial](https://jackbibby1.github.io/SCPA/articles/comparing_two_populations.html)
    for more detailed two group comparison with a specific data set
-   A
    [tutorial](https://jackbibby1.github.io/SCPA/articles/seurat_comparison.html)
    on how to use SCPA directly within a Seurat object
-   A
    [tutorial](https://jackbibby1.github.io/SCPA/articles/pseudotime.html)
    for multisample SCPA, comparing pathways across a pseudotime
    trajectory
-   A
    [tutorial](https://jackbibby1.github.io/SCPA/articles/systematic_tissue_comparison.html)
    of a systematic analysis of many cell types across multiple tissues
-   A
    [tutorial](https://jackbibby1.github.io/SCPA/articles/disease_comparison.html)
    for a systems level analysis of many cells types in disease
    (COVID-19)

## 3. Issues

To report any issues or submit comments please use:
<https://github.com/jackbibby1/SCPA/issues>
