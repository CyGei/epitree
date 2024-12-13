
# Check if a tree has a simple form:
# - one single introduction
# - all nodes are connected to at least one other node
check_tree <- function(tree) {
  # Check for a single introduction
  num_na <- sum(is.na(tree$from))
  if (num_na > 1) {
    stop("Invalid tree: More than one introduction (NAs in 'from').")
  }
  # Check for connectivity
  g <- graph_from_data_frame(d = tree, directed = TRUE)
  if (!is_connected(g, mode = "weak")) {
    stop("Invalid tree: The graph is not fully connected.")
  }
  return(TRUE)
}