
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCPA

<!-- badges: start -->
<!-- badges: end -->

SCPA is a method for pathway analysis in single cell RNA-seq data. It’s
based on a novel apprach to pathway analysis that defines pathway
activity as a change in multivariate distribution of a given pathway
across conditions, rather than enrichment or overrepresentation of
genes. This approach allows for a number of benefits over current
methods, but mainly: 1. You will identify pathways that show enrichment
in a given population AND also identify pathways with no overall
enrichment but alterations in the multivariate distribution of that
pathway. You essentially get the best of both worlds, as pathways with
changes in multivariate distribution but no overall enrichment are still
interestingly different pathways. 1. You can compare multiple conditions
simultaneously e.g. compare across 3 time points, or across multiple
phases of a pseodotime trajectory. This means you can assess pathway
activity through multiple stimulation times, or across cell
differentiation

## Installation

You can install the development version of SCPA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jackbibby1/SCPA")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SCPA)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
