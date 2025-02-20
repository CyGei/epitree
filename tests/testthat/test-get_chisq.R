make_tree_df <- function(n_cases, R, stochastic = TRUE) {
  make_tree(n_cases, R = R, stochastic) |>
    igraph::as_long_data_frame()
}

test_that("get_chisq performs basic chi-square test correctly", {
  chainA <- replicate(10, make_tree_df(20, 2), simplify = FALSE)
  chainB <- replicate(10, make_tree_df(20, 4), simplify = FALSE)

  result <- suppressWarnings(get_chisq(chainA, chainB))

  expect_s3_class(result, "htest")
  expect_true(!is.null(result$p.value))
  expect_equal(result$method, "Pearson's Chi-squared test")
})

test_that("get_chisq handles invalid inputs", {
  # Test with invalid data structure
  invalid_chain <- list(data.frame(x = 1, y = 2, z = 3))

  expect_error(get_chisq(invalid_chain))
})



#test that it is sensitive to id reshuffling
test_that("get_chisq is sensitive to id reshuffling", {
  chainA <- lapply(1:50, function(i) {
    make_tree(20, R = 2, stochastic = TRUE)
  })

  chainB <- lapply(1:50, function(i) {
    df <- epitree:::shuffle_graph_ids(chainA[[i]]) |>
      igraph::as_long_data_frame()
      subset(df, select = c("from", "to"))
  })

  chainA <- lapply(chainA, igraph::as_long_data_frame)

  result <- suppressWarnings(get_chisq(chainA, chainB))

  expect_true(result$p.value < 0.05) # Should be significantly different
})

#using very small replicates
test_that("get_chisq works with Fisher's test", {
  chainA <- replicate(4, make_tree_df(8, 2), simplify = FALSE)
  chainB <- replicate(4, make_tree_df(8, 4), simplify = FALSE)

  result <- get_chisq(chainA, chainB, method = "fisher")

  expect_s3_class(result, "htest")
  expect_match(result$method, "Fisher's Exact Test")
})

