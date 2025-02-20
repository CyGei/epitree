
library(testthat)
library(epitree)
library(igraph)

make_tree_df <- function(n_cases, R, stochastic = TRUE) {
  make_tree(n_cases, R = R, stochastic) |>
    igraph::as_long_data_frame()
}

test_that("compare_chains handles basic comparison correctly", {
  set.seed(123)
  chainA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)
  chainB <- replicate(10, make_tree_df(20, 4), simplify = FALSE)

  result <- suppressWarnings(compare_chains(chainA, chainB))

  expect_type(result, "list")
  expect_true(result$`Pr(>F)`[1] < 0.05)  # Should be significantly different
})

test_that("compare_chains works with identical chains", {
  set.seed(123)
  chainA <- replicate(10, make_tree_df(20, 3), simplify = FALSE)
  chainB <- replicate(10, make_tree_df(20, 3), simplify = FALSE)

  result <- suppressWarnings(compare_chains(chainA, chainB))

  expect_true(result$`Pr(>F)`[1] > 0.05)  # Should not be significantly different
})

test_that("compare_chains returns NA for f_tree = f_null", {
  chainA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)
  chainB <- replicate(10, make_tree_df(20, 4), simplify = FALSE)

  f_null <- function(tree) {
    return(rep(0, nrow(tree)))
  }

  result <- suppressWarnings(compare_chains(chainA, chainB, f_tree = f_null))
  expect_type(result, "list")
  expect_true(is.na(result$`Pr(>F)`[1]))  # Should be NA
})

test_that("compare_chains handles multiple chains", {
  chains <- lapply(1:4, function(i) {
    replicate(10, make_tree_df(20, 2), simplify = FALSE)
  })

  result <- suppressWarnings(do.call(compare_chains, chains))

  expect_type(result, "list")
  expect_true(result$`Pr(>F)`[1] > 0.05)  # Should not be significantly different
})

test_that("compare_chains handles invalid inputs appropriately", {
  chainA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)

  # Test with empty chain
  expect_error(compare_chains(list(), chainA))

  # Test with non-data frame elements
  invalid_chain <- list(1, 2, 3)
  expect_error(compare_chains(invalid_chain, chainA))

  # Test with data frames missing required columns
  invalid_df <- data.frame(x = 1:10, y = 1:10)
  invalid_chain <- replicate(10, invalid_df, simplify = FALSE)
  expect_error(compare_chains(invalid_chain, chainA))
})

test_that("compare_chains respects adonis2 arguments", {
  chainA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)
  chainB <- replicate(10, make_tree_df(20, 4), simplify = FALSE)

  result <- suppressWarnings(compare_chains(chainA, chainB,
                           adonis2_args = list(permutations = 10L)))

  permutation_call <- capture.output(print(result))[3]
  expect_match(permutation_call, "Number of permutations: 10")
})

#test that it is not sensitive to id reshuffling
test_that("compare_chains is not sensitive to id reshuffling", {
  set.seed(123)
  chainA <- lapply(1:50, function(i) {
    make_tree(20, R = 2, stochastic = TRUE)
  })

  chainB <- lapply(1:50, function(i) {
    df <- epitree:::shuffle_graph_ids(chainA[[i]]) |>
      igraph::as_long_data_frame()
    subset(df, select = c("from", "to"))
  })

  chainA <- lapply(chainA, igraph::as_long_data_frame)

  result <- suppressWarnings(compare_chains(chainA, chainB))

  expect_true(result$`Pr(>F)`[1] == 1)  # Should be exactly 1
})
