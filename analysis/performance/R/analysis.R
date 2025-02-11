# Set-Up ------------------------------------------------------------------
pacman::p_load(outbreaker2,
               ape,
               igraph,
               vegan,
               distcrete,
               tidyverse,
               scales,
               epitrix,
               furrr)
pacman::p_load_gh("CyGei/simulacr")
pacman::p_load_gh("CyGei/o2ools") # helper functions for outbreaker2
devtools::load_all()

source(here::here("analysis/performance/R/", "helpers.R"))

#run_date <- as.character(Sys.Date())
run_date <- "2025-02-09"
path <- "analysis/performance/data/"
dir.create(path = paste0(path, run_date),
           showWarnings = FALSE)


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
metadata <- data.frame(
  epidemic_size = epidemic_sizes[rep(seq_along(epidemic_sizes), each = n_simulations)],
  pair_IDX = seq_along(flat_pairs)
)


# Outbreak reconstruction -------------------------------------------------
flat_pairs <- flat_pairs[c(1,2, 101,102, 301, 302)]

future::plan("future::multisession", workers = 6)#future::availableCores() - 2
outs <- future_map(flat_pairs, function(pair) {
  A <- run_outbreaker(pair$A, ctd_fraction = 0.5)
  B <- run_outbreaker(pair$B, ctd_fraction = 0.65)

  list(A = A, B = A)
}, .options = furrr_options(seed = 123))
saveRDS(outs, file = paste0(path, run_date, "/outs.rds"))


# Process trees -----------------------------------------------------------
trees <- future_map(outs, function(out) {
  A <- out$A
  B <- out$B

  #Burnin
  A <- A[A$step > 1000, ]
  B <- B[B$step > 1000, ]

  tree_A <- o2ools::get_trees(A)
  tree_B <- o2ools::get_trees(B)

  tree_A <- lapply(tree_A, process_tree)
  tree_B <- lapply(tree_B, process_tree)

  return(list(A = tree_A, B = tree_B))
}, .options = furrr_options(seed = 123))


# Sensitivity Analysis ----------------------------------------------------
system.time({
  results <- future_map(trees, function(chain) {
    A <- chain$A
    B <- chain$B

    get_pvalues(
      chainA = A,
      chainB = B,
      sample_sizes = sample_sizes,
      overlap_freqs = overlap_freqs,
      n_repeats = 50
    )
  }, .options = furrr_options(seed = 123))
})

##################################
# Plot results
##################################
library(tidyverse, ggplot2, ggh4x)
df <- lapply(seq_along(results), function(i) {
  as.data.frame.table(results[[i]]) %>%
    mutate(epidemic_size = epidemic_sizes[i])
}) %>%
  bind_rows() %>%
  mutate(epidemic_size = factor(epidemic_size, levels = epidemic_sizes)) %>%
  rename(p_value = Freq)

alpha <- 0.05
annotations_df <- df %>%
  group_by(epidemic_size, sample_size, overlap_freq) %>%
  summarise(reject = mean(p_value < alpha) * 100,
            accept = mean(p_value >= alpha) * 100,
            .groups = "drop") %>%
  pivot_longer(
    cols = c("reject", "accept"),
    names_to = "category",
    values_to = "percentage"
  ) %>%
  mutate(
    row_title = "Epidemic Size",
    col_title = "Posterior Sample Size")

str(annotations_df)


p_props <- annotations_df %>%
  ggplot() +
  aes(x = percentage,
      y = factor(overlap_freq),
      fill = category) +
  ggh4x::facet_nested(
    rows = vars(row_title, epidemic_size),
    cols = vars(col_title, sample_size),
    scales = "fixed",
    space = "fixed",
    remove_labels = "none",
    nest_line = element_line(colour = "grey")
  )+
  geom_bar(stat = "identity",
           position = "stack",
           width = 0.7) +
  geom_text(
    aes(
      x = ifelse(category == "reject", percentage / 2, 100 - percentage / 2),
      label = ifelse(percentage > 0, round(percentage, 1), "")
    ),
    position = position_stack(vjust = 0.5),
    colour = "black",
    size = 3
  ) +
  scale_fill_manual(
    values =
      c(
        "reject" = "#FAD2CF",
        "accept" = "#CEEAD6"
      ),
    breaks = c("reject", "accept"),
    name = "P-value",
    labels = c("reject" = "Reject Null", "accept" = "Accept Null")
  ) +
  coord_cartesian(xlim = c(0, 100)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(family = "Fira Code",
                              #face = "bold",
                              size = 10),
    plot.title = element_text(size=12)
  ) +
  labs(x = "Percentage",
       y = "Overlap Frequency",
       title = "Proportion of p-values </> 0.05 by Overlap Frequency, Sample Size, and Epidemic Size")
p_props

#############################
# Test results
#############################
library(broom)
model_df <- df %>%
  mutate(
    category = ifelse(p_value < alpha, "reject", "accept"),
    reject = ifelse(p_value < alpha, 1, 0),
    accept = ifelse(p_value >= alpha, 1, 0),
  ) %>%
  select(epidemic_size,
         sample_size,
         overlap_freq,
         category,
         reject,
         accept)
#############################
# Chi-square test
#############################

chisq.test(table(model_df$category, model_df$epidemic_size))
chisq.test(table(model_df$category, model_df$sample_size))

#############################
# Logistic regression
#############################
# The logistic regression model is used to determine the effect of
# the overlap frequency and sample size on the probability of rejecting the null hypothesis.


glm(
  reject ~ sample_size,
  data = model_df,
  family = binomial(link = "logit")
) %>%
  broom::tidy() %>%
  mutate(OR = exp(estimate))

m1 <- glm(
  reject ~ sample_size + epidemic_size + as.numeric(as.character(overlap_freq)),
  data = model_df,
  family = binomial(link = "logit")
) %>%
  broom::tidy() %>%
  mutate(
    OR = exp(estimate),
    lower_ci = exp(estimate - 1.96 * std.error),
    upper_ci = exp(estimate + 1.96 * std.error)
  ) %>%
  filter(term != "(Intercept)")

ggplot(m1, aes(y = term, x = OR)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci),
                 height = 0.2) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "red") +
  scale_x_log10() +  # Log scale for better visualization of ORs
  labs(title = "Odds Ratios with 95% Confidence Intervals", x = "Odds Ratio (log scale)", y = "Variable") +
  theme_minimal()

# larger sample sizes strongly increase the odds of rejection
# epidemic size doesn't have an effect on the test outcome


# sensitivity model
model_df %>%
  filter(overlap_freq == 1) %>%
  glm(
    reject ~ sample_size + epidemic_size,
    data = .,
    family = binomial(link = "logit")
  ) %>%
  broom::tidy() %>%
  mutate(OR = exp(estimate))
