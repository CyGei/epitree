run_analysis_v2 <- function(trees_list,
                            sample_sizes = c(10, 50, 100, 500),
                            overlap_freqs = seq(0, 1, by = 0.2),
                            n_repeats = 100,
                            n_cores = future::availableCores() - 2,
                            seed = 123) {
  n_epidemic_sizes <- length(trees_list)

  # Modified get_pval to handle one epidemic size
  get_pval <- function(params,
                       chainA,
                       chainB,
                       n_repeats) {
    n <- params$sample_size
    pA <- params$overlap_freq
    pB <- 1 - pA

    lenA <- length(chainA)
    lenB <- length(chainB)
    total_chain <- c(chainA, chainB)

    # Probabilities
    probs <- c(rep(pA / lenA, lenA), rep(pB / lenB, lenB))

    # Pre-sample all at once
    mixed_samples <- replicate(n_repeats,
                               sample(
                                 x = total_chain,
                                 size = n,
                                 prob = probs,
                                 replace = TRUE
                               ),
                               simplify = FALSE)

    reference_samples <- replicate(n_repeats,
                                   sample(
                                     x = chainA,
                                     size = n,
                                     replace = TRUE
                                   ),
                                   simplify = FALSE)

    # Vectorize the test computation
    p_values <- mapply(function(ref, mix) {
      suppressWarnings(epitree::compare_chains(ref, mix)[, "Pr(>F)"][1])
    }, reference_samples, mixed_samples)

    return(p_values)
  }

  # Parameter grid
  grid <- expand.grid(
    sample_size = sample_sizes,
    overlap_freq = overlap_freqs,
    epidemic_size = seq_len(n_epidemic_sizes)
  )

  future::plan("future::multisession", workers = n_cores)

  p_values <- furrr::future_map(
    .x = split(grid, 1:nrow(grid)),
    .f = function(params) {
      i <- params$epidemic_size
      chainA <- trees_list[[i]][[1]]
      chainB <- trees_list[[i]][[2]]
      get_pval(params, chainA, chainB, n_repeats)
    },
    .options = furrr::furrr_options(
      seed = seed,
      globals = list(trees_list = trees_list)
    )
  )

  future::plan("future::sequential")

  # Reshape results into 4D array
  results <- array(
    data = unlist(p_values),
    dim = c(
      n_repeats,
      length(sample_sizes),
      length(overlap_freqs),
      n_epidemic_sizes
    ),
    dimnames = list(
      n_repeat = 1:n_repeats,
      sample_size = sample_sizes,
      overlap_freq = overlap_freqs,
      epidemic_size = seq_len(n_epidemic_sizes)
    )
  )

  return(results)
}
