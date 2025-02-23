
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixtree

<!-- badges: start -->

[![R-CMD-check](https://github.com/CyGei/mixtree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CyGei/mixtree/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `mixtree` package provides a statistical framework for comparing
sets of trees.

## Installation

You can install the development version of mixtree from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("CyGei/mixtree")
```

# Quick start

``` r
library(mixtree)
#> Registered S3 methods overwritten by 'adegraphics':
#>   method         from
#>   biplot.dudi    ade4
#>   kplot.foucart  ade4
#>   kplot.mcoa     ade4
#>   kplot.mfa      ade4
#>   kplot.pta      ade4
#>   kplot.sepan    ade4
#>   kplot.statis   ade4
#>   scatter.coa    ade4
#>   scatter.dudi   ade4
#>   scatter.nipals ade4
#>   scatter.pco    ade4
#>   score.acm      ade4
#>   score.mix      ade4
#>   score.pca      ade4
#>   screeplot.dudi ade4

# Simulate two chains of transmission trees
chainA <- lapply(1:100, function(i) {
  make_tree(20, R = 2, stochastic = TRUE) |>
    igraph::as_long_data_frame()
})
chainB <- lapply(1:100, function(i) {
  make_tree(20, R = 4, stochastic = TRUE) |>
    igraph::as_long_data_frame()
})

# Compare the two chains
result <- tree_test(chainA, chainB)
#> Warning in att$heading[2] <- deparse(match.call(), width.cutoff = 500L): number
#> of items to replace is not a multiple of replacement length
print(result)
#> Permutation test for adonis under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> (function (formula, data, permutations = 999, method = "bray", sqrt.dist = FALSE, add = FALSE, by = NULL, parallel = getOption("mc.cores"), na.action = na.fail, strata = NULL, ...) 
#>           Df SumOfSqs      R2      F Pr(>F)    
#> Model      1     7014 0.12621 28.599  0.001 ***
#> Residual 198    48558 0.87379                  
#> Total    199    55572 1.00000                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
