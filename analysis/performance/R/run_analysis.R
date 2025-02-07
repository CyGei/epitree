#' Run Sensitivity Analysis on Two Chains
#' This function runs a sensitivity analysis on two posterior chains of transmission trees (`chainA`, `chainB`) to check how the test
#' is affected by the sample size and the proportion of common trees in both chains.
#' @param chainA A vector representing the first chain.
#' @param chainB A vector representing the second chain.
#' @param sample_sizes A sequence of sample sizes to test. Default is \code{c(10, 50, 100, 500)}.
#' @param overlap_freqs A sequence of overlap proportions to test. Default is \code{seq(0, 1, by = 0.2)}.
#' @param n_repeats The number of repetitions to calculate the p-value. Default is 100.
#' @param n_cores The number of cores to use for parallel processing. Default is \code{future::availableCores() - 2}.
#' @param seed A seed for reproducibility. Default is 123.
#' @return A data frame with sample sizes, overlap indexes, and their corresponding p-values.
#' @keywords internal

run_analysis <- function(chainA,
                         chainB,
                         sample_sizes = c(10, 50, 100, 500),
                         overlap_freqs = seq(0, 1, by = 0.2),
                         n_repeats = 100,
                         n_cores = future::availableCores() - 2,
                         seed = 123) {
  # Parameter grid
  grid <- expand.grid(sample_size = sample_sizes, overlap_freq = overlap_freqs)

  get_pval <- function(params, chainA, chainB, n_repeats) {
    n <- params$sample_size
    pA <- params$overlap_freq
    pB <- 1 - pA

    lenA <- length(chainA)
    lenB <- length(chainB)
    total_chain <- c(chainA, chainB)

    # Probabilities
    probs <- c(rep(pA / lenA, lenA), rep(pB / lenB, lenB))

    # Pre-sample all at once
    mixed_samples <- replicate(n_repeats, sample(
      x = total_chain,
      size = n,
      prob = probs,
      replace = TRUE
    ), simplify = FALSE)

    reference_samples <- replicate(n_repeats, sample(
      x = chainA,
      size = n,
      replace = TRUE
    ), simplify = FALSE)

    # Vectorize the test computation
    p_values <- mapply(function(ref, mix) {
      suppressWarnings(epitree::compare_chains(ref, mix)[, "Pr(>F)"][1])
    }, reference_samples, mixed_samples)

    return(p_values)
  }

  future::plan("future::multisession", workers = n_cores)
  p_values <- furrr::future_map(
    .x = split(grid, 1:nrow(grid)),
    .f = ~ get_pval(.x, chainA, chainB, n_repeats),
    .options = furrr::furrr_options(
      seed = seed,
      globals = list(chainA = chainA, chainB = chainB)
    )
  )
  future::plan("future::sequential")

  results <- array(
    data = unlist(p_values),
    dim = c(n_repeats, length(sample_sizes), length(overlap_freqs)),
    dimnames = list(
      n_repeat = 1:n_repeats,
      sample_size = sample_sizes,
      overlap_freq = overlap_freqs
    )
  )

  return(results)
}
