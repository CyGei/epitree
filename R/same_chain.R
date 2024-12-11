# the purpose of this script is to compare two posterior sets of transmission trees sampled from the same chain.
pacman::p_load(outbreaker2, ape, igraph, vegan, distcrete)
pacman::p_load_gh("CyGei/o2ools")
source(here::here("R", "internals.R"))

linelist <- data.frame(
  ids = as.character(1:30),
  dates = outbreaker2::fake_outbreak$sample
)

set.seed(1)
chain <- outbreaker2::outbreaker(
  data = outbreaker2::outbreaker_data(
    dna = outbreaker2::fake_outbreak$dna,
    dates = outbreaker2::fake_outbreak$sample,
    ctd = outbreaker2::fake_outbreak$ctd,
    w_dens = outbreaker2::fake_outbreak$w,
    ids = linelist$ids
  ),
  # no unobserved cases should result in fully connected tree
  # without unobserved cases
  config = outbreaker2::create_config(
    init_pi = 1,
    move_pi = FALSE,
    find_import = FALSE
  )
)
chain <- chain[chain$step > 500, ]
plot(chain)

# Build trees ----------------------------------------------------
trees <- o2ools::get_trees(out = chain, ids = linelist$ids, kappa = FALSE)
trees <- lapply(trees, function(df) df[-1, ]) # remove introduction
invisible(lapply(trees, check_tree))


# Compare chains ----------------------------------------------------
source(here::here("R", "distances.R"))

set.seed(123)
chain1 <- sample(trees, size = 100, replace = TRUE)
chain2 <- sample(trees, size = 100, replace = TRUE)
compare_chains(chain1, chain2)
compare_chains(chain1, chain2, adonis2_args = list(permutations = 2000))

# Speed test ----------------------------------------------------
pacman::p_load(microbenchmark, profvis)
microbenchmark::microbenchmark(
  compare_chains(chain1, chain2, use_for = FALSE),
  compare_chains(chain1, chain2, use_for = TRUE)
)
# Unit: milliseconds
#                                             expr      min       lq     mean   median
#  compare_chains(chain1, chain2, use_for = FALSE) 523.1007 541.2627 563.1673 556.6668
#   compare_chains(chain1, chain2, use_for = TRUE) 488.7357 515.9245 542.6366 527.4246
#        uq      max neval
#  572.8273 737.4829   100
#  542.5292 827.6066   100

profvis::profvis({
  compare_chains(chain1, chain2)
})


# chisq test ----------------------------------------------------
pacman::p_load(tidyverse, magrittr)
ttabx <- sample(trees, size = 100, replace = TRUE) %$%
  bind_rows(.) %$%
  o2ools::ttable(from = .$from, to = .$to, levels = linelist$ids) 

ttaby <- sample(trees, size = 100, replace = TRUE) %$%
  bind_rows(.) %$%
  o2ools::ttable(from = .$from, to = .$to, levels = linelist$ids)

vegan::mantel(as.matrix(ttabx), as.matrix(ttaby))


get_chisq(
  chain_x = sample(trees, size = 100, replace = TRUE),
  chain_y = sample(trees, size = 100, replace = TRUE),
  levels = linelist$ids
)
