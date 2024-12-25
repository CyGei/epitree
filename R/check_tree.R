#' Check the Validity of a Transmission Tree
#'
#' Validates whether a given transmission tree meets the following topology:
#' - The graph is weakly connected (= a single introduction)
#' - The graph is a Directed Acyclic Graph (DAG) with no cycles or self-loops.
#' - Any node has at most one infector (no nodes have in-degree > 1).
#'
#' @param tree A data frame with two columns: `from` and `to`. The `from` column denotes the infector,
#'   and the `to` column denotes the infectee.
#'
#' @return Returns nothing if the tree meets all criteria. If the input is invalid,
#'   an error is thrown specifying the violated condition.
#'
#' @examples
#' # A valid tree
#' valid_tree <- data.frame(from = c("A", "A", "B", "C"), to = c("B", "C", "D", "E"))
#' check_tree(valid_tree)
#'
#' # Invalid tree: contains a cycle
#' invalid_tree_cycle <- data.frame(from = c("A", "B", "C"), to = c("B", "C", "A"))
#' \dontrun{
#' check_tree(invalid_tree_cycle)  # Error: must be a DAG
#' }
#'
#' # Invalid tree: multiple introductions
#' invalid_tree_multiple_roots <- data.frame(from = c("A", "B"), to = c("C", "D"))
#' \dontrun{
#' check_tree(invalid_tree_multiple_roots)  # Error: Must have exactly one root node
#' }
#'
#' @importFrom igraph graph_from_data_frame is_dag degree is_connected
#' @export

check_tree <- function(tree) {
  # Validate input
  if (!is.data.frame(tree) ||
      !all(c("from", "to") %in% colnames(tree))) {
    stop("Invalid tree: must be a data frame with columns 'from' and 'to'.")
  }

  if (sum(is.na(tree)) > 0) {
    stop("Invalid tree: cannot contain missing values.")
  }

  # Create directed graph
  g <- igraph::graph_from_data_frame(tree, directed = TRUE)

  # Perform checks
  if (!igraph::is_dag(g)) {
    stop("Invalid tree: must be a directed acyclic graph.")
  }

  if (!igraph::is_connected(g, mode = "weak")) {
    stop("Invalid tree: graph must be weakly connected.")
  }

  in_degrees <- igraph::degree(g, mode = "in")

  # if (sum(in_degrees == 0) != 1) {
  #   stop("Invalid tree: must have a single introduction.")
  # } # this is the same as weakly connected

  if (any(in_degrees > 1)) {
    stop("Invalid tree: each node must have at most one infector.")
  }


}
