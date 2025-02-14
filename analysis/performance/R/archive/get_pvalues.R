#' Run Sensitivity Analysis on Two Chains
#' This function runs a sensitivity analysis on two posterior chains of transmission trees (`chainA`, `chainB`) to check how the test
#' is affected by the sample size and the proportion of common trees in both chains.
#' @param chainA A list of transmission trees representing the first chain.
#' @param chainB A list of transmission trees representing the second chain.
#' @param sample_sizes A sequence of sample sizes to test. Default is \code{c(10, 50, 100, 500)}.
#' @param overlap_freqs A sequence of overlap proportions to test. Default is \code{seq(0, 1, by = 0.2)}.
#' @param n_repeats The number of repetitions to calculate the p-value. Default is 100.
#' @return A 3D array with the p-values for each combination of sample size and overlap proportion.
#' @keywords internal


get_pval <- function(chainA,
                     chainB,
                     overlap_freq,
                     sample_size,
                     n_repeats) {
  # Probabilities
  lenA <- length(chainA)
  lenB <- length(chainB)
  total_chain <- c(chainA, chainB)
  total_len <- lenA + lenB

  pA <- overlap_freq
  pB <- 1 - pA
  probs <- c(rep(pA / lenA, lenA), rep(pB / lenB, lenB))


  # Draw sample indices for all replicates in one go
  mix_idx <- matrix(
    sample.int(
      total_len,
      size = sample_size * n_repeats,
      prob = probs,
      replace = TRUE
    ),
    nrow = n_repeats
  )

  ref_idx <- matrix(sample.int(lenA, size = sample_size * n_repeats, replace = TRUE),
                    nrow = n_repeats)

  p_values <- numeric(n_repeats)
  for (i in seq_len(n_repeats)) {
    mixed_chain <- total_chain[mix_idx[i, ]]
    reference_chain <- chainA[ref_idx[i, ]]
    p_values[i] <- suppressWarnings(epitree::compare_chains(reference_chain, mixed_chain)[, "Pr(>F)"][1])
  }
  return(p_values)
}

get_pvalues <- function(chainA,
                        chainB,
                        sample_sizes = c(10, 50, 100, 500),
                        overlap_freqs = seq(0, 1, by = 0.2),
                        n_repeats = 50) {
  p_values <- array(
    NA,
    dim = c(n_repeats, length(sample_sizes), length(overlap_freqs)),
    dimnames = list(
      Replicate = 1:n_repeats,
      SampleSize = sample_sizes,
      OverlapFreq = overlap_freqs
    )
  )

  for (i in seq_along(sample_sizes)) {
    for (j in seq_along(overlap_freqs)) {
      p_values[, i, j] <- get_pval(
        chainA = chainA,
        chainB = chainB,
        overlap_freq = overlap_freqs[j],
        sample_size = sample_sizes[i],
        n_repeats = n_repeats
      )
    }
  }

  return(p_values)
}

