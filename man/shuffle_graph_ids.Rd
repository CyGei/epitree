% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_tree.R
\name{shuffle_graph_ids}
\alias{shuffle_graph_ids}
\title{Shuffle Node IDs in a Graph}
\usage{
shuffle_graph_ids(g, plot = FALSE)
}
\arguments{
\item{g}{An igraph object representing the graph.}

\item{plot}{Logical. If \code{TRUE}, the function will plot the shuffled graph. Default is \code{FALSE}.}
}
\value{
An igraph object with shuffled node IDs.
}
\description{
Randomly shuffles the IDs of the nodes in a given graph and optionally plots the shuffled graph.
}
\examples{
# Create an example graph
g <- make_tree(n_cases = 10, R = 2)

# Shuffle the node IDs
shuffled_graph <- shuffle_graph_ids(g, plot = TRUE)
}
