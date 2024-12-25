#' Perform Chi-Square or Fisher's test on posterior ancestries
#'
#' Test whether the distribution of infector-infectee pairs differs between chains.
#' @param ... Lists of data frames, where each data frame represents a transmission tree, with columns `from` (infector) and `to` (infectee). Two or more chains can be provided.
#' @param method A character string specifying the test to use: `"chisq"` (default) for Chi-Square
#'   or `"fisher"` for Fisher's Exact Test.
#' @param test_args Additional arguments passed to \code{\link[stats]{chisq.test}} or \code{\link[stats]{fisher.test}}. Must be a list.
#' @return A list of class `htest` containing the test results.
#'
#' @importFrom stats xtabs chisq.test fisher.test
#' @export

get_chisq <- function(...,
                      method = c("chisq", "fisher"),
                      test_args = list()) {
  method <- match.arg(method)

  make_tab <- function(chain) {
    df <- do.call(rbind, chain)
    if (ncol(df) != 2) {
      stop("Data must contain exactly two columns: 'from' and 'to'")
    }
    tab <- as.data.frame(table(df[, 1], df[, 2]))
    names(tab) <- c("from", "to", "Freq")
    return(tab)
  }
  chains <- list(...)
  tabs <- lapply(chains, make_tab)

  # Merge contingency tables
  tab <- Reduce(function(x, y)
    merge(x, y, by = c("from", "to"), all = TRUE), tabs)
  freq_cols <- grep("Freq", colnames(tab))
  tab <- tab[rowSums(tab[, freq_cols]) > 0, ]

  pairs <- paste(tab$from, tab$to, sep = "-")
  tab <- as.matrix(tab[, freq_cols])
  rownames(tab) <- pairs

  test_fn <- if (method == "fisher")
    stats::fisher.test
  else
    stats::chisq.test

  test_result <- do.call(test_fn, c(list(x = tab), test_args))

  # To print in regular format
  test_result$data.name <- "count data"

  return(test_result)
}

