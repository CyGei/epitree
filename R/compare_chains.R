#' Compare chains of posterior trees
#'
#' Performs a PERMANOVA (adonis2) to test for significant differences between posterior chains of transmission trees.
#'
#' @param ... Two or more posterior chains of transmission trees. Each chain is a list containing one or more data frames,
#' where each data frame represents a transmission tree with columns `from` (infector) and `to` (infectee).
#' @param f_tree A function to compute the pairwise distances between nodes within a tree.
#' Must take a single tree (data frame) as input and return a square matrix where dimensions match the number of nodes.
#' Default is `patristic`.
#' @param f_chain A function to compute the distance between two chains of transmission trees.
#' It must take two square matrices (produced by `f_tree`) as input and return a single numeric distance value.
#' Default is `euclidean`.
#' @param adonis2_args A list of additional arguments to pass to `vegan::adonis2`. Default is an empty list.
#'
#' @return A `vegan::adonis2` object containing the results of the PERMANOVA.
#'
#' @importFrom igraph graph_from_data_frame distances
#' @importFrom vegan adonis2
#' @importFrom stats as.dist
#' @importFrom utils combn
#'
#' @examples
#' set.seed(1)
#' #No difference in the chains
#' chainA <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 2, stochastic = TRUE)),
#'   simplify = FALSE)
#' chainB <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 2, stochastic = TRUE)),
#'   simplify = FALSE)
#' compare_chains(chainA, chainB)
#'
#' #Difference in the chains
#' chainC <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 4, stochastic = TRUE)),
#'   simplify = FALSE)
#' compare_chains(chainA, chainB, chainC)
#' @export

compare_chains <- function(...,
                           f_tree = patristic,
                           f_chain = euclidean,
                           adonis2_args = list()) {
  chains <- list(...)

  # Check for valid input
  if (length(chains) < 2) {
    stop("At least two chains are required.")
  }
  if (is.null(names(chains))) {
    names(chains) <- as.character(seq_along(chains))
  }

  # Generate chain_id for each data frame in the chains
  chain_lengths <- sapply(chains, function(chain) {
    ifelse(is.data.frame(chain), 1, length(chain))
  })
  chain_id <- factor(rep(names(chains), chain_lengths))

  if (all(vapply(chains, is.data.frame, logical(1)))) {
    master_chain <- chains # if only 1 dataframe per chain
  } else {
    master_chain <- unlist(chains, recursive = FALSE)
  }

  # Validate each tree in master_chain (may remove)
  for (i in seq_along(master_chain)) {
    validate_tree(master_chain[[i]])
  }

  within_list <- lapply(master_chain, f_tree)

  # Get all unique pair combinations#
  n <- length(master_chain)
  pairs <- t(utils::combn(n, 2))
  d <- matrix(0, n, n)


  distances <- apply(pairs, 1, function(p) {
    f_chain(within_list[[p[1]]], within_list[[p[2]]])
  })
  d[lower.tri(d)] <- distances

  # adonis2
  if (!is.list(adonis2_args)) {
    stop("adonis2_args must be a list.")
  }

  adonis <- do.call(
    vegan::adonis2,
    c(
      list(
        formula = stats::as.dist(d) ~ chain_id,
        data = data.frame(chain_id = chain_id)
      ),
      adonis2_args
    )
  )

  # Warning in att$heading[2] <- deparse(match.call(), width.cutoff = 500L) :
  # number of items to replace is not a multiple of replacement length
  return(adonis)
}
