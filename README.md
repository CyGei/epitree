
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
print(result)
#> Permutation test for adonis under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> (function (formula, data, permutations = 999, method = "bray", sqrt.dist = FALSE, add = FALSE, by = NULL, parallel = getOption("mc.cores"), na.action = na.fail, strata = NULL, ...) 
#>           Df SumOfSqs      R2      F Pr(>F)    
#> Model      1     6730 0.13248 30.236  0.001 ***
#> Residual 198    44071 0.86752                  
#> Total    199    50801 1.00000                  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
