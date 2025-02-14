pacman::p_load(tidyverse, ggh4x)
run_date <- "2025-02-09"
path <- "analysis/performance/data/"
# Load data ---------------------------------------------------------------
metadata <- left_join(x = readRDS(paste0(path, run_date, "/grid.rds")),
                  y = readRDS(paste0(path, run_date, "/sim_metadata.rds")),
                  by = "tree_id")

p_values <- readRDS(paste0(path, run_date, "/p_values.rds"))

results <- tibble::tibble(metadata, p_value = p_values) %>%
  mutate(overlap_freq = as.factor(overlap_freq),
         epidemic_size = as.factor(epidemic_size),
         sample_size = as.factor(sample_size)) %>%
  tidyr::unnest_longer(p_value,indices_to = "replicate")

# plot --------------------------------------------------------------------
alpha <- 0.05
annotations_df <- results %>%
  group_by(epidemic_size, sample_size, overlap_freq) %>%
  summarise(
    reject = mean(p_value < alpha) * 100,
    accept = mean(p_value >= alpha) * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c("reject", "accept"),
    names_to = "category",
    values_to = "percentage"
  ) %>%
  mutate(row_title = "Epidemic Size", col_title = "Posterior Sample Size")

str(annotations_df)


p_props <- annotations_df %>%
  ggplot() +
  aes(x = overlap_freq,
      y = percentage,
      fill = category) +
  ggh4x::facet_nested(
    rows = vars(row_title, epidemic_size),
    cols = vars(col_title, sample_size),
    scales = "fixed",
    space = "fixed",
    remove_labels = "none",
    nest_line = element_line(colour = "grey")
  ) +
  geom_bar(stat = "identity",
           position = "stack",
           width = 0.7) +
  geom_text(
    aes(
      y = ifelse(category == "reject", percentage / 2, 100 - percentage / 2),
      label = ifelse(percentage > 0, round(percentage, 1), "")
    ),
    position = position_stack(vjust = 0.5),
    colour = "black",
    size = 3
  ) +
  scale_fill_manual(
    values =
      c("reject" = "#FAD2CF", "accept" = "#CEEAD6"),
    breaks = c("reject", "accept"),
    name = "P-value",
    labels = c("reject" = "Reject Null", "accept" = "Accept Null")
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(family = "Fira Code", # face = "bold",
                              size = 10),
    plot.title = element_text(size = 12)
  ) +
  labs(x = "Overlap Frequency", y = "Percentage", title = "Proportion of p-values </> 0.05 by Overlap Frequency, Sample Size, and Epidemic Size")
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


glm(reject ~ sample_size,
    data = model_df,
    family = binomial(link = "logit")) %>%
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
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "red") +
  scale_x_log10() + # Log scale for better visualization of ORs
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
