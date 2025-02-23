library(testthat)
library(mixtree)
library(igraph)

make_tree_df <- function(n_cases, R, stochastic = TRUE) {
  make_tree(n_cases, R = R, stochastic) |>
    igraph::as_long_data_frame()
}

#####################
# PERMANOVA CHECKS
#####################
test_that("tree_test correctly detects differences using PERMANOVA", {
  set.seed(123)
  setA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)
  setB <- replicate(10, make_tree_df(20, 4), simplify = FALSE)

  result <- suppressWarnings(tree_test(setA, setB))

  expect_type(result, "list")
  # Expect a significant difference since the sets have different R values
  expect_true(result$`Pr(>F)`[1] < 0.05)
})

test_that("tree_test returns non-significant result when comparing identical sets", {
  set.seed(123)
  setA <- replicate(10, make_tree_df(20, 3), simplify = FALSE)
  setB <- replicate(10, make_tree_df(20, 3), simplify = FALSE)

  result <- suppressWarnings(tree_test(setA, setB))

  # Expect no significant difference between identical sets
  expect_true(result$`Pr(>F)`[1] > 0.05)
})

test_that("tree_test returns NA when within_dist function returns no variability", {
  setA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)
  setB <- replicate(10, make_tree_df(20, 4), simplify = FALSE)

  f_null <- function(tree) {
    return(rep(0, nrow(tree)))
  }

  result <- suppressWarnings(tree_test(setA, setB, within_dist = f_null))
  expect_type(result, "list")
  # When within_dist provides no variability, the test statistic should be NA
  expect_true(is.na(result$`Pr(>F)`[1]))
})

test_that("tree_test correctly handles more than two sets", {
  sets <- lapply(1:4, function(i) {
    replicate(10, make_tree_df(20, 2), simplify = FALSE)
  })

  result <- suppressWarnings(do.call(tree_test, sets))

  expect_type(result, "list")
  # When comparing multiple identical sets, no significant difference is expected
  expect_true(result$`Pr(>F)`[1] > 0.05)
})

test_that("tree_test errors for invalid input formats", {
  setA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)

  # Expect error when an empty set is provided
  expect_error(tree_test(list(), setA))

  # Expect error when the set contains elements that are not data frames
  invalid_set <- list(1, 2, 3)
  expect_error(tree_test(invalid_set, setA))

  # Expect error when data frames lack the required 'from' and 'to' columns
  invalid_df <- data.frame(x = 1:10, y = 1:10)
  invalid_set <- replicate(10, invalid_df, simplify = FALSE)
  expect_error(tree_test(invalid_set, setA))
})

test_that("tree_test correctly passes adonis2 arguments to vegan::adonis2", {
  setA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)
  setB <- replicate(10, make_tree_df(20, 4), simplify = FALSE)

  result <- suppressWarnings(tree_test(setA, setB,
                                       test_args = list(permutations = 10L)))

  permutation_call <- capture.output(print(result))[3]
  # Check that the number of permutations specified is correctly used
  expect_match(permutation_call, "Number of permutations: 10")
})

# Test that the function is robust to changes in tree ID ordering
test_that("tree_test is robust to reshuffling of tree IDs", {
  set.seed(123)
  setA <- lapply(1:50, function(i) {
    make_tree(20, R = 2, stochastic = TRUE)
  })

  setB <- lapply(1:50, function(i) {
    df <- mixtree:::shuffle_graph_ids(setA[[i]]) |>
      igraph::as_long_data_frame()
    subset(df, select = c("from", "to"))
  })

  setA <- lapply(setA, igraph::as_long_data_frame)

  result <- suppressWarnings(tree_test(setA, setB))

  # When IDs are reshuffled, the test should indicate no difference (p-value exactly 1)
  expect_true(result$`Pr(>F)`[1] == 1)
})

#####################
# CHISQ TEST
#####################

test_that("tree_test performs the basic chi-square test correctly", {
  setA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)
  setB <- replicate(10, make_tree_df(20, 4), simplify = FALSE)

  result <- suppressWarnings(tree_test(setA, setB, method = "chisq"))

  expect_s3_class(result, "htest")
  expect_true(!is.null(result$p.value))
  # Check that the method label indicates a Pearson's Chi-squared test
  expect_equal(result$method, "Pearson's Chi-squared test")
})

# Test that chi-square method detects differences when IDs are reshuffled
test_that("tree_test detects differences in chi-square test with ID reshuffling", {
  setA <- lapply(1:50, function(i) {
    make_tree(20, R = 2, stochastic = TRUE)
  })

  setB <- lapply(1:50, function(i) {
    df <- mixtree:::shuffle_graph_ids(setA[[i]]) |>
      igraph::as_long_data_frame()
    subset(df, select = c("from", "to"))
  })

  setA <- lapply(setA, igraph::as_long_data_frame)

  result <- suppressWarnings(tree_test(setA, setB, method = "chisq"))

  # Expect a significant difference due to ID reshuffling affecting pair frequencies
  expect_true(result$p.value < 0.05)
})

# Test using small replicate numbers with Fisher's Exact Test
test_that("tree_test performs Fisher's Exact Test correctly with small replicates", {
  setA <- replicate(4, make_tree_df(8, 2), simplify = FALSE)
  setB <- replicate(4, make_tree_df(8, 4), simplify = FALSE)

  result <- tree_test(setA, setB, method = "fisher")

  expect_s3_class(result, "htest")
  # Confirm that the method label corresponds to Fisher's Exact Test
  expect_match(result$method, "Fisher's Exact Test")
})

