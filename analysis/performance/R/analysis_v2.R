##################################
# Set-up
##################################
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
source(here::here("analysis/performance/R/run_analysis.R"))
run_date <- as.character(Sys.Date())
run_date <- "2025-02-03"
path <- "analysis/performance/data/"
dir.create(path = paste0(path, run_date),
           showWarnings = FALSE)

##################################
# Helpers
##################################
# Run a simulation until the desired epidemic size is reached
simulate_until <- function(epidemic_size, R_values, duration = 100) {
  population_size <- round(epidemic_size / epitrix::R02AR(mean(R_values)))
  repeat {
    sim <- simulacr::simulate_outbreak(
      duration = duration,
      population_size = population_size,
      R_values = R_values,
      dist_incubation = outbreaker2::fake_outbreak$w,
      dist_generation_time = outbreaker2::fake_outbreak$w
    )$data %>%
      relabel_tree(id_cols = c("id", "source"), date_col = "date_onset")

    if (nrow(sim) == epidemic_size)
      return(sim)
  }
}
# Add a day for cases with the same onset date as the first case
shift_init_date <- function(sim,
                            date_col = "date_onset",
                            shift = 1) {
  dates <- sim[[date_col]]
  min_date <- min(dates, na.rm = TRUE)
  idx <- which(dates == min_date)
  # Shift all but the first occurrence
  if (length(idx) > 1) {
    sim[[date_col]][idx[-1]] <- sim[[date_col]][idx[-1]] + shift
  }
  return(sim)
}

# Relabel the trees based on the date of onset
label_ids <- function(id, date) {
  ordered_ids <- id[order(date)]
  new_labels <- as.character(seq_along(ordered_ids))
  setNames(new_labels, ordered_ids)[id]
}

relabel_tree <- function(df, id_cols, date_col = "date_onset") {
  ids <- label_ids(df[[id_cols[1]]], df[[date_col]]) # Use the first column for labelling

  # Apply the relabelling all specified id columns
  for (col in id_cols) {
    df[[col]] <- ids[df[[col]]]
  }

  return(df)
}



##################################
# Grid parameters
##################################
sample_sizes <- c(20, 50, 100, 200)
overlap_freqs <- seq(0, 1, by = 0.2)
epidemic_sizes <- c(20, 50, 100, 200)




##################################
# Simulations by epidemic size
##################################

future::plan("future::multisession", workers = length(epidemic_sizes))
sims <- future_map(epidemic_sizes, function(epidemic_size) {
  sim_ss    <- simulate_until(epidemic_size, rnbinom(n = 100, size = 0.2, mu = 3))
  sim_linear <- simulate_until(epidemic_size, 1L)
  list(sim_ss = sim_ss, sim_linear = sim_linear)
}, .options = furrr_options(seed = 123))
saveRDS(sims, file = paste0(path, run_date, "/sims.rds"))

##################################
# Outbreak reconstruction
##################################
outs <- furrr::future_map(seq_along(epidemic_sizes), function(i) {
  epidemic_size <- epidemic_sizes[i]

  # Extract the simulated outbreaks for this epidemic size.
  # Ensure that the first case has a unique onset date
  sim_ss     <- shift_init_date(sims[[i]]$sim_ss)
  sim_linear <- shift_init_date(sims[[i]]$sim_linear)

  # Add contact data for improved accuracy
  ctd_sample1 <- sample(2:epidemic_size, round(epidemic_size * 0.5))
  ctd_sample2 <- sample(2:epidemic_size, round(epidemic_size * 0.65))

  config <- outbreaker2::create_config(
    init_tree = "star",
    move_mu = FALSE,
    init_pi = 1,
    move_pi = FALSE,
    init_kappa = 1,
    move_kappa = FALSE,
    find_import = FALSE)


  out_ss <- outbreaker2::outbreaker(
    data = outbreaker2::outbreaker_data(
      dates = sim_ss$date_onset,
      w_dens = outbreaker2::fake_outbreak$w,
      f_dens = outbreaker2::fake_outbreak$w,
      ctd = sim_ss[ctd_sample1, c("source", "id")],
      ids = sim_ss$id
    ),
    config = config
  )
  out_ss <- identify(out_ss, sim_ss$id)

  out_linear <- outbreaker2::outbreaker(
    data = outbreaker2::outbreaker_data(
      dates = sim_linear$date_onset,
      w_dens = outbreaker2::fake_outbreak$w,
      f_dens = outbreaker2::fake_outbreak$w,
      ctd = sim_linear[ctd_sample2, c("source", "id")],
      ids = sim_linear$id
    ),
    config = config
  )
  out_linear <- identify(out_linear, sim_linear$id)

  list(out_ss = out_ss, out_linear = out_linear)
}, .options = furrr::furrr_options(seed = 123))
saveRDS(outs, file = paste0(path, run_date, "/outs.rds"))


##################################
# Process trees
##################################
trees <- lapply(seq_along(epidemic_sizes), function(i) {
  epidemic_size <- epidemic_sizes[i]
  out_ss <- outs[[i]]$out_ss
  out_linear <- outs[[i]]$out_linear

  #Burnin
  out_ss <- out_ss[out_ss$step > 1000, ]
  out_linear <- out_linear[out_linear$step > 1000, ]


  # get the trees
  tree_ss <- o2ools::get_trees(out_ss, sims[[i]]$sim_ss$id)
  tree_linear <- o2ools::get_trees(out_linear, sims[[i]]$sim_linear$id)

  # Process the trees
  tree_ss <- lapply(tree_ss, process_tree)
  tree_linear <- lapply(tree_linear, process_tree)

  list(tree_ss = tree_ss, tree_linear = tree_linear)
})


##################################
# Run analysis for each epidemic size
##################################
system.time({
  results <- lapply(seq_along(epidemic_sizes), function(i) {
    epidemic_size <- epidemic_sizes[i]
    trees_ss <- trees[[i]]$tree_ss
    trees_linear <- trees[[i]]$tree_linear

    result <- run_analysis(
      chainA = trees_ss,
      chainB = trees_linear,
      sample_sizes = sample_sizes,
      overlap_freqs = overlap_freqs,
      n_repeats = 100
    )
    return(result)
  })
})

saveRDS(results, file = paste0(path, run_date, "/results.rds"))
#70790.74 sec ~ 20hrs



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
