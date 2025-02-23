#' Perform PERMANOVA on Sets of Transmission Trees
#'
#' Tests for significant differences between sets of transmission trees using PERMANOVA (via \code{vegan::adonis2}).
#'
#' @param ... Two or more sets of transmission trees. Each set is a list of dataframes with columns \code{from} (infector) and \code{to} (infectee).
#' @param within_dist A function to compute pairwise distances within a tree. Takes a dataframe, returns a square matrix. Default is \code{\link{patristic}}.
#' @param between_dist A function to compute distance between two trees. Takes two matrices, returns a numeric value. Default is \code{\link{euclidean}}.
#' @param test_args A list of additional arguments to pass to \code{vegan::adonis2}. Default is an empty list.
#'
#' @return A \code{vegan::adonis2} object containing the test results.
#'
#' @importFrom igraph graph_from_data_frame distances
#' @importFrom vegan adonis2
#' @importFrom stats as.dist
#' @importFrom utils combn
#'
#' @examples
#' set.seed(1)
#' # No difference in the sets
#' setA <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 2, stochastic = TRUE)
#' ),
#' simplify = FALSE
#' )
#' setB <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 2, stochastic = TRUE)
#' ),
#' simplify = FALSE
#' )
#' mixtree:::permanova_test(setA, setB)
#'
#' # Difference in the sets
#' setC <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 4, stochastic = TRUE)
#' ),
#' simplify = FALSE
#' )
#' mixtree:::permanova_test(setA, setB, setC)
#' @keywords internal

permanova_test <- function(...,
                           within_dist = patristic,
                           between_dist = euclidean,
                           test_args = list()) {
  sets <- list(...)
  validate_sets(sets)

  # Tag each tree with a set ID
  set_lengths <- sapply(sets, length)
  set_id <- factor(rep(seq_along(sets), set_lengths))
  master_set <- unlist(sets, recursive = FALSE)

  # Validate all trees (@NOTE: may remove)
  for (i in seq_along(master_set)) {
    validate_tree(master_set[[i]])
  }

  # Compute within-tree distance
  within_list <- lapply(master_set, within_dist)

  # Get all unique pair combinations
  n <- length(master_set)
  pairs <- t(utils::combn(n, 2))
  d <- matrix(0, n, n)

  # Compute between-tree distance
  distances <- apply(pairs, 1, function(p) {
    between_dist(within_list[[p[1]]], within_list[[p[2]]])
  })
  d[lower.tri(d)] <- distances

  # PERMANOVA
  if (!is.list(test_args)) {
    stop("test_args must be a list.")
  }

  adonis <- do.call(
    vegan::adonis2,
    c(
      list(
        formula = stats::as.dist(d) ~ set_id,
        data = data.frame(set_id = set_id)
      ),
      test_args
    )
  )

  # Warning in att$heading[2] <- deparse(match.call(), width.cutoff = 500L) :
  # number of items to replace is not a multiple of replacement length
  return(adonis)
}
