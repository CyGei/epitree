#' Validate a Transmission Tree
#'
#' Checks if a transmission tree meets specific topology criteria for our test.
#' The tree must be a directed acyclic graph (DAG), weakly connected, and have at most one infector per node.
#'
#' @param tree A data frame with columns `from` and `to` representing the transmission tree.
#'
#' @return Invisible `TRUE` if the tree is valid. Throws an error if invalid.
#'
#' @importFrom igraph graph_from_data_frame is_dag is_connected degree
#'
#' @examples
#' good_tree <- data.frame(from = c(1, 2, 3), to = c(2, 3, 4))
#' validate_tree(good_tree)
#' bad_tree <- data.frame(from = c(1, 2, 3), to = c(2, 3, 2))
#' try(validate_tree(bad_tree))
#' @export
validate_tree <- function(tree) {
  # Validate input
  if (!is.data.frame(tree) ||
    !all(c("from", "to") %in% colnames(tree))) {
    stop("must be a data frame with columns 'from' and 'to'.")
  }
  if (any(is.na(tree))) {
    stop("cannot contain missing values.")
  }
  g <- igraph::graph_from_data_frame(tree, directed = TRUE)

  # Perform checks
  if (!igraph::is_dag(g)) {
    stop("must be a directed acyclic graph.")
  }
  if (!igraph::is_connected(g, mode = "weak")) {
    stop("graph must be weakly connected.") # i.e. all nodes reachable from root
  }
  if (any(igraph::degree(g, mode = "in") > 1)) {
    stop("each node must have at most one infector.")
  }

  invisible(TRUE)
}


#' Validate a Set of Transmission Trees
#'
#' Ensures that the input is a list containing at least one dataframe.
#'
#' @param set A list containing at least one dataframe.
#' @return Invisible `TRUE` if the set is valid. Throws an error if invalid.
#' @keywords internal

validate_set <- function(set) {
  # Check if it's strictly a list (not a data frame)
  if (!inherits(set, "list")) {
    stop("set must be a list of dataframes, not a single dataframe or other type.")
  }

  # Check if it's empty
  if (length(set) == 0) {
    stop("set is an empty list.")
  }

  # Check all elements are dataframes
  is_df <- vapply(set, is.data.frame, logical(1))
  if (!all(is_df)) {
    first_non_df <- which(!is_df)[1] # Find first offending element
    stop(sprintf("Element %d in chain must be a data frame.", first_non_df))
  }

  invisible(TRUE)
}

#' Validate sets of transmission trees
#'
#' Checks that the provided input is a list of at least two valid sets of transmission trees.
#' Each set is expected to be a list containing at least one data frame, as verified by
#' \code{\link{validate_set}}.
#'
#' @param sets A list where each element represents a set of transmission trees. Each set must be a list
#' containing one or more data frames.
#'
#' @return Invisible `TRUE` if the sets are valid. Throws an error if invalid.
#'
#' @details
#' At least two sets are provided.
#' Each set is a list (and not a data frame itself).
#' Each set contains at least one element.
#' Every element in each set is a data frame.
#'
#' @seealso \code{\link{validate_set}} for validating an individual set.
#' @keywords internal
validate_sets <- function(sets) {
  if (length(sets) < 2) {
    stop("At least two chains are required.")
  }
  for (i in seq_along(sets)) {
    validate_set(sets[[i]])
  }

  invisible(TRUE)
}

#' Remove the Introduction from a Transmission Tree
#'
#' Will remove the single introduction (row with an `NA` in the `from` column).
#'
#' @param tree A data frame with columns `from` and `to` representing the transmission tree.
#' @return A cleaned tree with the introduction removed.
#' @noRd
remove_intro <- function(tree) {
  total_nas <- sum(is.na(tree))
  from_nas <- sum(is.na(tree$from))

  # Remove the introduction
  if (total_nas == 1 && from_nas == 1) {
    tree <- tree[!is.na(tree$from), ]
  } else {
    stop("invalid tree: must contain exactly one introduction.")
  }
}

#' Process a Transmission Tree
#'
#' Will attempt to validate the tree.
#' If the tree is invalid, it will remove the introduction and try again.
#'
#' @param tree A data frame with columns `from` and `to`.
#' @return The validated and potentially cleaned tree.
#' @noRd
process_tree <- function(tree, silent = FALSE) {
  # Attempt to validate the tree first
  tryCatch(
    {
      validate_tree(tree)
      return(tree)
    },
    error = function(e) {
      tree <- remove_intro(tree)
      validate_tree(tree)
      return(tree)
    }
  )
}
