#' Generate a Transmission Tree
#'
#' Creates a transmission tree with a specified number of cases and branches per case.
#' The tree can be generated with fixed or Poisson-distributed branching factors.
#'
#' @param n_cases Integer. The total number of cases (nodes) in the tree.
#' @param R Integer. The fixed number of branches per case when \code{stochastic} is \code{FALSE},
#'           or the mean of the Poisson distribution when \code{stochastic} is \code{TRUE}.
#' @param stochastic Logical. If \code{TRUE}, the number of branches per case is sampled from
#'                   a Poisson distribution with mean \code{R}. Default is \code{FALSE}.
#' @param plot Logical. If \code{TRUE}, the function will plot the generated tree. Default is \code{FALSE}.
#' @return An igraph object representing the transmission tree.
#' @importFrom igraph make_graph layout_as_tree V as_long_data_frame
#' @importFrom stats rpois
#' @examples
#' # Generate a deterministic transmission tree
#' deterministic_tree <- make_tree(n_cases = 15, R = 2, stochastic = FALSE, plot = TRUE)
#'
#' # Generate a stochastic transmission tree
#' random_tree <- make_tree(n_cases = 15, R = 2, stochastic = TRUE, plot = TRUE)
#' @export

make_tree <- function(n_cases, R = 2, stochastic = FALSE, plot = FALSE) {
  edges <- c()
  current_nodes <- 1
  next_node <- 2

  while(next_node <= n_cases) {
    new_nodes <- c()
    for(node in current_nodes) {
      # Determine the branching factor
      if(stochastic) {
        R_current <- stats::rpois(1, lambda = R)
      } else {
        R_current <- R
      }

      for(i in 1:R_current) {
        if(next_node > n_cases) {
          break
        }
        # Add an edge from the current node to the next node
        edges <- c(edges, node, next_node)
        new_nodes <- c(new_nodes, next_node)
        next_node <- next_node + 1
      }
    }
    current_nodes <- new_nodes
    if(length(current_nodes) == 0) {
      break
    }
  }

  # Create the graph object
  g <- igraph::make_graph(edges, directed = TRUE)

  validate_tree(igraph::as_long_data_frame(g))

  # Plot the graph if requested
  if(plot) {
    plot(
      g,
      vertex.label = igraph::V(g)$name,
      layout = igraph::layout_as_tree,
      vertex.color = ifelse(stochastic, "lightgreen", "skyblue"),
      edge.arrow.size = 0.5,
      main = ifelse(stochastic, "Randomised Transmission Tree", "Deterministic Transmission Tree")
    )
  }

  return(g)
}



#' Shuffle Node IDs in a Graph
#'
#' Randomly shuffles the IDs of the nodes in a given graph and optionally plots the shuffled graph.
#'
#' @param g An igraph object representing the graph.
#' @param plot Logical. If \code{TRUE}, the function will plot the shuffled graph. Default is \code{FALSE}.
#' @return An igraph object with shuffled node IDs.
#' @importFrom igraph vcount permute V layout_as_tree
#' @examples
#' # Create an example graph
#' g <- make_tree(n_cases = 10, R = 2)
#'
#' # Shuffle the node IDs
#' shuffled_graph <- shuffle_graph_ids(g, plot = TRUE)
#' @export
shuffle_graph_ids <- function(g, plot = FALSE) {
  # Get the number of vertices in the graph
  num_vertices <- igraph::vcount(g)

  # Generate a random permutation of vertex IDs
  permutation <- sample(1:num_vertices, num_vertices)

  # Apply the permutation to the graph
  g_shuffled <- igraph::permute(g, permutation)

  # Update vertex labels to reflect the new IDs
  V(g_shuffled)$name <- as.character(1:num_vertices)

  if(plot){
    plot(
      g_shuffled,
      vertex.label = igraph::V(g_shuffled)$name,
      layout = igraph::layout_as_tree,
      vertex.color = "orange",
      edge.arrow.size = 0.5,
      main = "Shuffled Transmission Tree"
    )
  }

  return(g_shuffled)
}
