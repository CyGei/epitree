library(testthat)
library(igraph)

test_that("check_tree function works correctly", {
  # Valid tree
  valid_tree <- data.frame(from = c("A", "A", "B", "C"), to = c("B", "C", "D", "E"))
  expect_silent(check_tree(valid_tree))

  # Invalid tree: contains a cycle
  invalid_tree_cycle <- data.frame(from = c("A", "B", "C"), to = c("B", "C", "A"))
  expect_error(check_tree(invalid_tree_cycle), "must be a directed acyclic graph")

  # Invalid tree: not weakly connected (i.e. more than 1 introduction)
  invalid_tree_disconnected <- data.frame(from = c("A", "B", "D"), to = c("B", "C", "E"))
  expect_error(check_tree(invalid_tree_disconnected), "graph must be weakly connected")

  # Invalid tree: node has multiple infectors
  invalid_tree_multiple_infectors <- data.frame(from = c("A", "B", "B", "A"), to = c("B", "C", "D", "D"))
  expect_error(check_tree(invalid_tree_multiple_infectors), "each node must have at most one infector")

})

