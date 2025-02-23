#' Test Differences Between Sets of Transmission Trees
#'
#' Performs a statistical test to assess whether there are significant differences between sets of transmission trees.
#' Supports PERMANOVA (via \code{"vegan::adonis2"}), Chi-Square, or Fisher's Exact Test.
#'
#' @param ... Two or more sets of transmission trees. Each set must be a list of data frames with columns \code{from} (infector) and \code{to} (infectee).
#' @param method A character string specifying the test method. Options are \code{"permanova"}, #' \code{"chisq"}, or \code{"fisher"}. Default is \code{"permanova"}.
#' @param within_dist A function to compute pairwise distances within a tree for PERMANOVA. Takes a data frame, returns a square matrix. Default is \code{\link{patristic}}.
#' @param between_dist A function to compute distance between two trees for PERMANOVA. Takes two matrices, returns a numeric value. Default is \code{\link{euclidean}}.
#' @param test_args A list of additional arguments to pass to the underlying test function (\code{vegan::adonis2}, \code{stats::chisq.test}, or \code{stats::fisher.test}). Default is an empty list.

#' @return
#' - For \code{"permanova"}: A \code{"vegan::adonis2"} object containing the test results.
#' - For \code{"chisq"} or \code{"fisher"}: An \code{"htest"} object with the test results.
#'
#' @details
#' This function compares sets of transmission trees using one of three statistical tests.
#'
#' **PERMANOVA**: Evaluates whether the topological distribution of transmission trees differs between sets.
#'   - **Null Hypothesis (H0)**: Transmission trees in all sets are drawn from the same distribution, implying similar topologies.
#'   - **Alternative Hypothesis (H1)**: At least one set of transmission trees comes from a different distribution.
#'
#' **Chi-Square or Fisher’s Exact Test**: Evaluates whether the distribution of infector-infectee pairs differs between sets.
#'   - **Null Hypothesis (H0)**: The frequency of infector-infectee pairs is consistent across all sets.
#'   - **Alternative Hypothesis (H1)**: The frequency of infector-infectee pairs differs between at least two sets.
#'
#' @importFrom igraph graph_from_data_frame distances
#' @importFrom vegan adonis2
#' @importFrom stats as.dist chisq.test fisher.test complete.cases
#' @importFrom utils combn
#'
#' @examples
#' set.seed(1)
#' # Generate example sets
#' setA <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 2, stochastic = TRUE)
#' ), simplify = FALSE)
#' setB <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 2, stochastic = TRUE)
#' ), simplify = FALSE)
#' setC <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 4, stochastic = TRUE)
#' ), simplify = FALSE)
#'
#' # PERMANOVA test
#' tree_test(setA, setB, setC,  method = "permanova")
#'
#' # Chi-Square test
#' tree_test(setA, setB, setC, method = "chisq")
#' @export

tree_test <- function(
    ...,
    method = c("permanova", "chisq", "fisher"),
    within_dist = patristic,
    between_dist = euclidean,
    test_args = list()) {

  method <- match.arg(method)

  if (method == "permanova") {
    result <- permanova_test(...,
      within_dist = within_dist,
      between_dist = between_dist,
      test_args = test_args
    )
  } else {
    result <- chisq_test(...,
      method = method,
      test_args = test_args
    )
  }

  return(result)
}
