
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


# computes the adjacency matrix
# indicate whether pairs of nodes are adjacent or not in the graph
adj_matrix <- function(from, to, levels = NULL, method = "adjacency", diag = FALSE, ...) {
  if (!is.vector(from) || !is.vector(to) || length(from) !=
    length(to)) {
    stop("'from' and 'to' must be vectors of the same length.")
  }
  from <- as.character(from)
  to <- as.character(to)
  if (is.null(levels)) {
    levels <- unique(sort(c(from, to)))
  } else {
    if (!is.vector(levels) || any(is.na(levels))) {
      stop("'levels' must be a vector of non-missing values.")
    }
    if (!all(unique(sort(c(from, to))) %in% levels)) {
      stop("The unique values in 'from' and 'to' must be a subset of 'levels'.")
    }
    levels <- unique(levels)
  }
  if (length(levels) < 2) {
    stop("There must be at least two group levels in the data")
  }

  if (method == "adjacency") {
    m <- table(c(from, to), c(to, from)) > 0 # in case of cycles - avoid overcounting
    m <- 1 * m
  }
  # if (method == "directional"){
  # ttable is not applicable because it accounts for directionality (from -> to)
  #   m <- o2ools::ttable(from, to, levels = levels)
  # }
  d <- as.dist(m, ...)

  return(d)
}

# # Check if all ttables have the same number of rows
# same_row_count <- all(sapply(ttables, nrow) == nrow(ttables[[1]]))
# # Check if all ttables have the same order of from-to pairs
# same_order <- all(sapply(ttables, function(df) {
#   identical(df[, c("from", "to")], ttables[[1]][, c("from", "to")])
# }))





# Draft code ----------------------------------------------------

# for more efficient computation
#t(combn(length(master_chain), 2))

# iGraphs 
# graphs <- lapply(trees, \(x) igraph::graph_from_data_frame(x, directed = TRUE))

# visNetwork::visIgraph(graphs[[1]]) |>
#   visNetwork::visIgraphLayout(layout = "layout_nicely") |>
#   visNetwork::visOptions(nodesIdSelection = TRUE, highlightNearest = TRUE)
# dists <- lapply(graphs, \(x) distances(x, mode = "all", algorithm = "unweighted"))

# epicontacts
# epic <- make_epicontacts(
#   linelist = linelist, contacts = trees[[1]],
#   directed = TRUE
# )
# summary(epic)
# plot(epic)

# distances(
#   g,
#   mode = "out"
# )