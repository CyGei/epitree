#' Run Sensitivity Analysis on Two Chains
#' This function runs a sensitivity analysis on two posterior chains of transmission trees (`chainA`, `chainB`) to check how the test
#' is affected by the sample size and the proportion of common trees in both chains.
#' @param chainA A vector representing the first chain.
#' @param chainB A vector representing the second chain.
#' @param sample_sizes A sequence of sample sizes to test. Default is \code{seq(10, 1000, by = 200)}.
#' @param overlap_freqs A sequence of overlap proportions to test. Default is \code{seq(0, 1, by = 0.2)}.
#' @param n_repeats The number of repetitions to calculate the p-value. Default is 100.
#' @param n_cores The number of cores to use for parallel processing. Default is \code{future::availableCores() - 2}.
#' @param seed A seed for reproducibility. Default is 123.
#' @return A data frame with sample sizes, overlap indexes, and their corresponding p-values.
#' @keywords internal

run_analysis <- function(chainA,
                         chainB,
                         sample_sizes = seq(10, 1000, by = 200),
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

    # Probabilities
    lenA <- length(chainA)
    lenB <- length(chainB)
    probs <- rep(c(pA, pB), times = c(lenA, lenB))
    probs <- probs / sum(probs)

    replicate(n_repeats, {
      # Generate samples
      mixed <- sample(
        x = c(chainA, chainB),
        size = n,
        prob = probs,
        replace = TRUE
      )

      reference <- sample(x = chainA,
                          size = n,
                          replace = TRUE)

      # Return p-value
      suppressWarnings(epitree::compare_chains(reference, mixed)[, "Pr(>F)"][1])
    })
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
