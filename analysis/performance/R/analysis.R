#################################################################
run_date <- as.character(Sys.Date())
dir.create(path = paste0("analysis/data/", run_date),
           showWarnings = FALSE)
system.time({
  # simulation for type 2 error
  result <- run_analysis(
    chainA = chain_uni,
    chainB = chain_ss,
    sample_sizes = c(50, 200, 1000),
    overlap_freqs = seq(0, 1, 0.1),
    #generate_sequence(0, 1, 0.1)
    n_repeats = 1000
  )
  saveRDS(result, file = paste0("analysis/data/", run_date, "/result.rds"))
}) -> time
saveRDS(time, file = paste0("analysis/data/", run_date, "/time.rds"))

#################################################################

library(tidyverse)
# date_folder <- list.files("analysis/data", full.names = TRUE) |>
#   vapply(\(x) as.numeric(as.Date(file.info(x)$ctime)), numeric(1)) |>
#   which.max() %>%
#   names()

result <-
  readRDS(here::here("analysis/data/", run_date, "/result.rds")) %>%
  reshape2::melt() %>%
  mutate(across(.cols = c(sample_size, overlap_freq), .fns = as.factor))


ggplot(result, aes(x = overlap_freq, , y = value, color = sample_size)) +
  facet_wrap( ~ sample_size) +
  geom_violin() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "bottom")

str(result)
length(levels(result$overlap_freq))))
result %>%
  mutate(overlap_freq = as.numeric(as.character(overlap_freq))) %>%
  ggplot(aes(
    x = sample_size,
    y = value,
    color = overlap_freq,
    group = interaction(sample_size, overlap_freq)
  )) +
  geom_violin() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "bottom")

pacman::p_load(gghalves)
result %>%
  mutate(title = "sample size") %>%
  ggplot(aes(
    y = value,
    x = overlap_freq,
    group = interaction(sample_size, overlap_freq, title),
  )) +
  facet_nested( ~ title + sample_size, scales = "free", space = "free") +  # Add hierarchical facets
  gghalves::geom_half_violin(
    scale = "width",
    fill = "grey",
    color = "black",
    alpha = 0.5,
    draw_quantiles = c(0.05, 1 - 0.05),
    side = "r",
    nudge = 0.025
  ) +
  gghalves::geom_half_boxplot(side = "l",
                              outlier.size = 0.1,
                              nudge = 0.025) +
  #stat_binline(bins = 50, scale = 0.9, draw_baseline = FALSE)+
  geom_hline(yintercept = 0.05,
             linetype = "dashed",
             colour = "orange") +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      size = 0.5
    ),
    strip.text.x = element_text(face = "bold")
  )+
  labs(
    x = "Overlap Frequency",
    y = "p-value",
    title = "Distribution of p-values by Sample Size and Overlap Frequency"
  )



# Calculate confusion matrix metrics for each combination
metrics <- result %>%
  mutate(
    true_label = case_when(overlap_freq == 0 ~ 1, overlap_freq == 1 ~ 0, TRUE ~ NA_real_),
    pred_label = ifelse(value <= 0.05, 1, 0)
  ) %>%
  filter(overlap_freq == 0 | overlap_freq == 1) %>%
  group_by(sample_size, overlap_freq) %>%
  summarise(
    TP = sum(pred_label == 1 & true_label == 1),
    TN = sum(pred_label == 0 & true_label == 0),
    FP = sum(pred_label == 1 & true_label == 0),
    # Type I error
    FN = sum(pred_label == 0 & true_label == 1),
    # Type II error

    # Sensitivity (True Positive Rate)
    sensitivity = TP / sum(true_label == 1),

    # Specificity (True Negative Rate)
    specificity = TN / sum(true_label == 0),
  )

# ROC points
