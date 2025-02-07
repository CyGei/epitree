install.packages("epichains")
library(epichains)
set.seed(32)
# Simulate chains
k <- 0.2 # Dispersion parameter
R <- 2 # Basic reproduction number
prob <- k / (k + R)

sim_chains <- simulate_chains(
  n_chains = 1,
  statistic = "size",
  offspring_dist = rnbinom,
  size = k,
  prob = prob,
  stat_threshold = 50L,
  pop = 100L
)
as.data.frame(sim_chains)

longest_chain <- sim_chains[sim_chains$chain == which.max(
  unname(table(sim_chains$chain))
), ]

# Convert to data.frame to view the whole data
as.data.frame(longest_chain)

make_epichains <- function(chain){
  epicontacts::make_epicontacts(
    linelist = subset(chain, select = -infector),
    contacts = subset(chain, select = c(infectee, infector)),
    na_rm_contacts = TRUE,
    directed = TRUE
  )
}

plot_epichains <- function(epicontacts, x_axis = "time"){
  plot(epicontacts, x_axis = x_axis)
}

sim_chains %>% make_epichains() %>% plot_epichains()

make_linear_tree <- function(n_cases) {
  tree <- data.frame(
    id = as.character(seq_len(n_cases)),
    source = c(NA, as.character(seq_len(n_cases - 1))),
    R = c(rep(1, n_cases - 1), 0),
    date_onset = seq_len(n_cases)
  )

  return(tree)
}
