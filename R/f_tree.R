#' Compute the Patristic distance matrix
#'
#' The patristic distance is the number of generations separating any two nodes in a transmission tree.
#'
#' @param tree A data frame representing a transmission tree, with the first column containing the infector IDs and the second the infectee IDs.
#' @return A square, symmetric matrix of patristic distances between nodes.
#' @examples
#' tree <- data.frame(from = c(1, 1, 2, 2, 3, 3), to = c(2, 3, 4, 5, 6, 7))
#' patristic(tree)
#' @import igraph
#' @export
patristic <- function(tree) {
  g <- igraph::graph_from_data_frame(tree, directed = FALSE)
  igraph::distances(g)
}

#' Compute the Abouheif distance matrix
#'
#' The Abouheif distance is the product of the number of direct descendants of each node in the path between two nodes.
#' It is a measure of the number of transmission events between two nodes.
#'
#' @param tree A data frame representing a transmission tree, with the first column containing the infector IDs and the second the infectee IDs.
#' @return A square, symmetric matrix of Abouheif distances between nodes.
#' @examples
#' tree <- data.frame(from = c(1, 1, 2, 2, 3, 3), to = c(2, 3, 4, 5, 6, 7))
#' abouheif(tree)
#' @import igraph
#' @export
abouheif <- function(tree) {
  #@TODO: ask if we include node1 and node2 in P?

  g <- igraph::graph_from_data_frame(tree, directed = FALSE)
  node_ids <- igraph::V(g)$name
  node_n <- length(node_ids)
  abouheif_dist <- matrix(0, nrow = node_n, ncol = node_n)
  rownames(abouheif_dist) <- colnames(abouheif_dist) <- node_ids

  # identify pairs of nodes with > 1 patristic distance
  # i.e. where there's at least one node between them
  patristic_dist <- igraph::distances(g)
  pairs <- which(lower.tri(patristic_dist) &
                   patristic_dist > 1, arr.ind = TRUE)

  # degree of each node (DD = R_case + 2)
  DD <- igraph::degree(g)

  # identify the nodes in the path (P) between pairs
  DDP <- apply(pairs, 1, function(p) {
    nodes <- igraph::shortest_paths(graph = g,
                                    from = p[[1]],
                                    to = p[[2]])$vpath[[1]]
    nodes <- nodes[-c(1, length(nodes))]
    prod(DD[nodes])
  })

  abouheif_dist[pairs] <- DDP

  #return symmetric matrix
  abouheif_dist <- abouheif_dist + t(abouheif_dist)

  return(abouheif_dist)
}


#' Compute the Kendall distance matrix
#'
#' Kendall's distance measures the depth of the most recent common infector (MRCI) for each pair of nodes with respect to the source (patient 0).
#'
#' @param tree A data frame representing a transmission tree, with the first column containing the infector IDs and the second the infectee IDs.
#' @return A square, symmetric matrix of Kendall's distances between nodes.
#' @examples
#' tree <- data.frame(from = c(1, 1, 2, 2, 3, 3), to = c(2, 3, 4, 5, 6, 7))
#' kendall(tree)
#' @importFrom treespace findMRCIs
#' @references
#' A Metric to Compare Transmission Trees - M Kendall · 2018
#' @seealso \code{\link[treespace]{findMRCIs}}
#' @export
kendall <- function(tree) {
  treespace::findMRCIs(as.matrix(tree))$mrciDepths
}
