pacman::p_load_gh("reconhub/simulacr")



#################################################################
dir.create(path = paste0("analysis/data/", Sys.Date()),
           showWarnings = FALSE)
system.time({
  # simulation for type 2 error
  result <- run_analysis(
    chainA = trees[[1]],
    chainB = trees[[2]],
    sample_sizes = c(50, 200, 1000),
    overlap_indexes = seq(0, 1, 0.1) #generate_sequence(0, 1, 0.1)
  )
  saveRDS(result, file = paste0("analysis/data/", Sys.Date(), "/result.rds"))
}) -> time
saveRDS(time, file = paste0("analysis/data/", Sys.Date(), "/time.rds"))

#################################################################

library(tidyverse)
date_folder <- list.files("analysis/data", full.names = TRUE) |>
  vapply(\(x) as.numeric(as.Date(file.info(x)$ctime)), numeric(1)) |>
  which.max() %>%
  names()

results <- bind_rows(
  readRDS(paste0(date_folder, "/type1.rds")) %>%
    mutate(chains = "identical"),
  readRDS(paste0(date_folder, "/type2.rds")) %>%
    mutate(chains = "different")
) %>%
  mutate(across(
    .cols = c(sample_size, overlap_index),
    .fns = as.factor
  )) %>%
  mutate(
    true_label = ifelse(chains == "different", 1, 0),
    pred_label = ifelse(p_value <= 0.05, 1, 0)
  )


# Calculate confusion matrix metrics for each combination
metrics <- results %>%
  group_by(sample_size, overlap_index) %>%
  summarise(
    TP = sum(pred_label == 1 & true_label == 1),
    TN = sum(pred_label == 0 & true_label == 0),
    FP = sum(pred_label == 1 & true_label == 0), # Type I error
    FN = sum(pred_label == 0 & true_label == 1), # Type II error

    # Sensitivity (True Positive Rate)
    sensitivity = TP / sum(true_label == 1),

    # Specificity (True Negative Rate)
    specificity = TN / sum(true_label == 0),
  )

ggplot(metrics, aes(x = sample_size, y = overlap_index)) +
  geom_tile(aes(fill = sensitivity), colour = "white") +
  scale_fill_gradient(low = "grey", high = "#008080") +
  labs(
    title = "Sensitivity",
    x = "Sample Size",
    y = "Overlap Index"
  ) +
  theme_minimal()
