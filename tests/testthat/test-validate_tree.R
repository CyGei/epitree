library(testthat)
library(igraph)

# validate_tree
test_that("validate_tree works correctly with valid trees", {
  valid_tree <- data.frame(from = c("A", "A", "B", "C"), to = c("B", "C", "D", "E"))
  expect_silent(validate_tree(valid_tree))
})

test_that("validate_tree throws errors for invalid trees", {

  #Tree with a cycle
  cycle_tree <- data.frame(from = c("A", "B", "C"), to = c("B", "C", "A"))
  expect_error(validate_tree(cycle_tree), "must be a directed acyclic graph")

  #Tree with an introduction
  intro_tree <- data.frame(from = c(NA, "A", "B"), to = c("A", "B", "C"))
  expect_error(validate_tree(intro_tree), "cannot contain missing values")

  #Tree weakly connected
  weak_tree <- data.frame(from = c("A", "C", "D"), to = c("B", "D", "E"))
  expect_error(validate_tree(weak_tree), "graph must be weakly connected")

  # Tree with nodes having multiple infectors
  multi_infector_tree <- data.frame(from = c("A", "B", "B"), to = c("C", "C", "D"))
  expect_error(validate_tree(multi_infector_tree), "each node must have at most one infector")
})

#process_tree
test_that("process_tree works correctly with valid trees", {
  valid_tree <- data.frame(from = c("A", "A", "B", "C"), to = c("B", "C", "D", "E"))
  #Valid tree should be returned without modification
  result <- process_tree(valid_tree)
  expect_equal(result, valid_tree)
})

test_that("process_tree prepares and validates intially invalid trees", {
  # Remove introduction and validate
  intro_tree <- data.frame(from = c(NA, "A", "B"), to = c("A", "B", "C"))
  result <- process_tree(intro_tree)
  expected <- data.frame(from = c("A", "B"), to = c("B", "C"))
  rownames(result) <- NULL
  rownames(expected) <- NULL
  expect_equal(result, expected)

})
