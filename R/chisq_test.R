#' Perform Chi-Square Test on Sets of Transmission Trees
#'
#' Tests whether the distribution of infector-infectee pairs differs between sets of transmission trees.
#'
#' @param ... Two or more sets of transmission trees. Each set is a list of data frames with columns \code{from} and \code{to}.
#' @param method Test to use: \code{"chisq"} for Chi-Square or \code{"fisher"} for Fisher's Exact Test. Default is \code{"chisq"}.
#' @param test_args A list of additional arguments for \code{stats::chisq.test} or \code{stats::fisher.test}. Default is an empty list.
#'
#' @return An \code{htest} object with the test results.
#'
#' @importFrom stats chisq.test fisher.test complete.cases
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
#' mixtree:::chisq_test(setA, setB)
#'
#' # Difference in the sets
#' setC <- replicate(10, igraph::as_long_data_frame(
#'   make_tree(n_cases = 10, R = 4, stochastic = TRUE)
#' ),
#' simplify = FALSE
#' )
#' mixtree:::chisq_test(setA, setB, setC)
#' @keywords internal

chisq_test <- function(...,
                       method = c("chisq", "fisher"),
                       test_args = list()) {
  method <- match.arg(method)
  sets <- list(...)
  validate_sets(sets)

  make_tab <- function(set) {
    df <- do.call(rbind, set)
    if (ncol(df) != 2) {
      stop("Data must contain exactly two columns: 'from' and 'to'")
    }
    tab <- as.data.frame(table(df[, 1], df[, 2]), stringsAsFactors = FALSE)
    names(tab) <- c("from", "to", "Freq")
    return(tab)
  }
  tabs <- lapply(sets, make_tab)

  # Merge contingency tables
  tab <- Reduce(function(x, y) {
    merge(x, y, by = c("from", "to"), all = TRUE)
  }, tabs)
  freq_cols <- grep("Freq", colnames(tab))
  tab <- tab[stats::complete.cases(tab), ]
  tab <- tab[rowSums(tab[, freq_cols]) > 0, ]

  pairs <- paste(tab$from, tab$to, sep = "-")
  tab <- as.matrix(tab[, freq_cols])
  rownames(tab) <- pairs

  test_fn <- if (method == "fisher") {
    stats::fisher.test
  } else {
    stats::chisq.test
  }

  test_result <- do.call(test_fn, c(list(x = tab), test_args))

  # To print in regular format
  test_result$data.name <- "count data"

  return(test_result)
}
