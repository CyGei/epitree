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
#' @examples
#' good_tree <- data.frame(from = c(1,2,3), to = c(2,3,4))
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
    stop("graph must be weakly connected.") #i.e. all nodes reachable from root
  }
  if (any(igraph::degree(g, mode = "in") > 1)) {
    stop("each node must have at most one infector.")
  }

  invisible(TRUE)
}



#' Remove the Introduction from a Transmission Tree
#'
#' Will remove the single introduction (row with an `NA` in the `from` column).
#'
#' @param tree A data frame with columns `from` and `to` representing the transmission tree.
#'
#' @return A cleaned tree with the introduction removed.
#' @keywords internal
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
#'
#' @return The validated and potentially cleaned tree.
#' @keywords internal
process_tree <- function(tree, silent = FALSE) {
  # Attempt to validate the tree first
  tryCatch({
    validate_tree(tree)
    return(tree)
  }, error = function(e) {
    tree <- remove_intro(tree)
    validate_tree(tree)
    return(tree)
  })
}
