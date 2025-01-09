run_chisq_analysis <- function(chainA,
                               chainB,
                               sample_sizes = seq(10, 1000, by = 200),
                               overlap_freqs = seq(0, 1, by = 0.2),
                               n_repeats = 100,
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

  p_values <- lapply(split(grid, 1:nrow(grid)),
                     ~ get_pval(.x, chainA, chainB, n_repeats),)

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



result <- run_chisq_analysis(
  chainA = chain_uni,
  chainB = chain_ss,
  sample_sizes = c(50, 200, 1000),
  overlap_freqs = seq(0, 1, 0.1), #generate_sequence(0, 1, 0.1)
  n_repeats = 1000
)
