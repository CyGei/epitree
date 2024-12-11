# Compare two posterior sets of transmission trees sampled from 2 different chains.
pacman::p_load(outbreaker2, ape, igraph, vegan, distcrete, truncnorm)
pacman::p_load_gh("CyGei/o2ools")
source(here::here("R", "internals.R"))

linelist <- data.frame(
  ids = as.character(1:30),
  dates = outbreaker2::fake_outbreak$sample
)

gt <- lapply( list(
  distcrete("norm", interval = 1, mean = 30, sd = 1, w = 0.5),
  distcrete("norm", interval = 1, mean = 10, sd = 2, w = 0.5)
), function(d) d$d(1:50) )


set.seed(1)
chains <-
  lapply(1:2, function(i){
    x <- outbreaker2::outbreaker(
      data = outbreaker2::outbreaker_data(
        dna = outbreaker2::fake_outbreak$dna,
        dates = outbreaker2::fake_outbreak$sample,
        ctd = outbreaker2::fake_outbreak$ctd,
        w_dens = gt[[i]],
        ids = linelist$ids
      ),
      config = outbreaker2::create_config(
        init_pi = 1,
        move_pi = FALSE,
        find_import = FALSE
      )
    )
    x <- x[x$step > 500, ]
  })


# Build trees ----------------------------------------------------

trees <- lapply(chains, function(chain) {
  t <- o2ools::get_trees(out = chain, ids = linelist$ids, kappa = FALSE)
  t <- lapply(t, function(df) df[-1, ]) # remove introduction
  invisible(lapply(t, check_tree))
  return(t)
  })



# Compare chains ----------------------------------------------------
source(here::here("R", "distances.R"))

chain1 <- trees[[1]]
chain2 <- trees[[2]]

set.seed(123)
compare_chains(chain1, chain2)

# Evaluate how the proportion of shared trees affects the p-value by n_trees
n_trees <- c(10, 75, 150)
freq_common  <- seq(0, 1, by = 0.05)
results <- do.call(rbind, lapply(n_trees, function(n) {

  ref_chain <- trees[[1]][sample(1:length(trees[[1]]), n)]
  mixed_chain <- trees[[2]][sample(1:length(trees[[2]]), n)]

  # Compute p-values for each freq_common
  data.frame(
    n_trees = rep(n, length(freq_common)),
    freq_common = freq_common,
    p_value = sapply(freq_common, function(perc) {
      n_replace <- round(perc * n)

      # Create the mixed chain (that contains freq_common trees from the reference chain)
      new_mixed <- c(
        sample(mixed_chain, n - n_replace, replace = FALSE),
        sample(ref_chain, n_replace, replace = FALSE)
      )
      new_mixed <- sample(new_mixed, length(new_mixed))

      adonis_result <- compare_chains(ref_chain, new_mixed)
      adonis_result[1, "Pr(>F)"]
    })
  )
}))

pacman::p_load(ggplot2, scales)
ggplot(results,
       aes(
         x = freq_common,
         y = p_value,
         color = as.character(n_trees),
         group = as.character(n_trees)
       )) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "Proportion of trees from reference chain", y = "p-value") +
  scale_y_log10() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10),
                     labels = scales::percent_format()) +
  theme_bw()


