# Set-Up ------------------------------------------------------------------
pacman::p_load(
  outbreaker2,
  ape,
  igraph,
  vegan,
  distcrete,
  tidyverse,
  scales,
  epitrix,
  furrr
)
pacman::p_load_gh("CyGei/simulacr")
pacman::p_load_gh("CyGei/o2ools") # helper functions for outbreaker2
devtools::load_all()

source(here::here("analysis/performance/R/", "helpers.R"))
source(here::here("analysis/performance/R/", "get_pval.R"))
# run_date <- as.character(Sys.Date())
run_date <- "2025-02-09"
path <- "analysis/performance/data/"
dir.create(
  path = paste0(path, run_date),
  showWarnings = FALSE
)


# Grid Parameters ---------------------------------------------------------
sample_sizes <- c(20, 50, 100, 200)
overlap_freqs <- seq(0, 1, by = 0.2)
epidemic_sizes <- c(20, 50, 100, 200)


# Outbreak Simulations -----------------------------------------------------
# Parameters
R_values_A <- rnbinom(100, size = 0.1, mu = 2)
R_values_B <- rpois(100, lambda = 2)
tolerance <- 0.2
n_simulations <- 100
max_attempts <- n_simulations * 100
try_gain <- 0.7

future::plan("future::multisession", workers = length(epidemic_sizes))
sims <- future_map(epidemic_sizes, function(N) {
  A <- simulate_outbreaks(
    target_size = N,
    R_values = R_values_A,
    tolerance = tolerance,
    n_simulations = n_simulations,
    max_attempts = max_attempts,
    try_gain = try_gain
  )

  B <- simulate_outbreaks(
    target_size = N,
    R_values = R_values_B,
    tolerance = tolerance,
    n_simulations = n_simulations,
    max_attempts = max_attempts,
    try_gain = try_gain
  )
  return(list(A = A, B = B))
}, .options = furrr_options(seed = 123))
saveRDS(sims, file = paste0(path, run_date, "/sims.rds"))

# Pair simulations together
pairs <- future_map(seq_along(epidemic_sizes), function(i) {
  matching_pairs(sims[[i]][["A"]], sims[[i]][["B"]], n = n_simulations)
}, .options = furrr_options(seed = 123))
flat_pairs <- unlist(pairs, recursive = FALSE)

# Record metadata
sim_metadata <- data.frame(
  epidemic_size = epidemic_sizes[rep(seq_along(epidemic_sizes), each = n_simulations)],
  tree_id = as.character(seq_along(flat_pairs))
)
#saveRDS(sim_metadata, file = paste0(path, run_date, "/sim_metadata.rds"))

# Outbreak reconstruction -------------------------------------------------
future::plan("future::multisession", workers = 6) # future::availableCores() - 2
outs <- future_map(flat_pairs, function(pair) {
  A <- run_outbreaker(pair$A, ctd_fraction = 0.5)
  B <- run_outbreaker(pair$B, ctd_fraction = 0.65)

  list(A = A, B = A)
}, .options = furrr_options(seed = 123))
saveRDS(outs, file = paste0(path, run_date, "/outs.rds"))


# Process trees -----------------------------------------------------------
future::plan("future::multisession", workers = future::availableCores() - 2)
trees <- future_map(outs, function(out) {
  A <- out$A
  B <- out$B

  # Burnin
  A <- A[A$step > 1000, ]
  B <- B[B$step > 1000, ]

  tree_A <- o2ools::get_trees(A)
  tree_B <- o2ools::get_trees(B)

  tree_A <- lapply(tree_A, epitree:::process_tree)
  tree_B <- lapply(tree_B, epitree:::process_tree)

  return(list(A = tree_A, B = tree_B))
}, .options = furrr_options(seed = 123))
saveRDS(trees, file = paste0(path, run_date, "/trees.rds"))

# save each pair of tree individually in a trees folder with index as id
# trees <- readRDS(paste0(path, run_date, "/trees.rds"))
# for(i in seq_along(trees)){
#   file <- paste0(path, run_date,"/trees/tree_", i, ".rds")
#   saveRDS(trees[[i]], file = file)
# }


# Sensitivity Analysis ----------------------------------------------------
files <- list.files(paste0(path, run_date, "/trees"))
grid <- expand.grid(
  tree_id = as.character(seq_along(files)),
  overlap_freq = seq(0, 1, by = 0.2),
  sample_size = c(20, 50, 100, 200)
)
saveRDS(grid, file = paste0(path, run_date, "/grid.rds"))

future::plan("future::multisession", workers = future::availableCores() - 2)
system.time({
  p_values <- future_map(asplit(grid, 1), function(param) {
    tree <- readRDS(paste0(path, run_date, "/trees/tree_", param[["tree_id"]], ".rds"))
    get_pval(
      chainA = tree[["A"]],
      chainB = tree[["B"]],
      overlap_freq = as.double(param[["overlap_freq"]]),
      sample_size = as.integer(param[["sample_size"]]),
      n_repeats = 50
    )
  }, .options = furrr_options(seed = 123))
})
p_values <- saveRDS(p_values, file = paste0(path, run_date, "/p_values.rds"))
