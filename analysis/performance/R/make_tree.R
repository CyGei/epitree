make_superspreaders <- function(tree, ss_ids, R_lambdas) {
  if (!all(ss_ids %in% tree$id)) stop("One or more invalid case IDs")

  updated_tree <- tree

  # Get superspreaders' positions (indices)
  ss_positions <- match(ss_ids, updated_tree$id)

  # Generate the number of secondary infections, ensuring they don’t exceed available cases
  max_cases <- nrow(updated_tree) - ss_positions
  n_secondary_cases <- pmin(rpois(length(ss_ids), R_lambdas), max_cases)

  # Update R values for superspreaders
  updated_tree$R[ss_positions] <- n_secondary_cases

  # Assign secondary cases
  for (i in seq_along(ss_positions)) {
    ss_pos <- ss_positions[i]
    n_cases <- n_secondary_cases[i]

    if (n_cases == 0) next  # Skip if no secondary infections

    infectee_positions <- (ss_pos + 1):(ss_pos + n_cases)
    updated_tree$source[infectee_positions] <- ss_ids[i]

    # Maintain continuity by updating the next case's source
    next_pos <- ss_pos + n_cases + 1
    if (next_pos <= nrow(updated_tree)) {
      updated_tree$source[next_pos] <- updated_tree$id[ss_pos + n_cases]
    }
  }

  return(updated_tree)
}


#######################
# make_tree
#######################
# Function to generate a tree with a given number of cases and R value,
# The logic ensures that each case infects R other cases, if possible.
# As the number of cases is finite, the actual R value may be lower.

# n_cases: number of cases to generate
# R_value: average number of secondary cases per primary case
# Returns a data frame with columns id, source, R, and date_onset

make_tree <- function(n_cases, R_value) {
  # Initialise the tree
  tree <- data.frame(
    id = as.character(seq_len(n_cases)),
    source = rep(NA, n_cases),
    R = rep(0, n_cases)
  )

  # Track the next available case
  next_available <- 2

  # Force Case 1 to infect at least one person
  max_possible <- n_cases - next_available + 1
  R_1 <- min(max(1, rpois(1, R_value)), max_possible)  # Ensure at least 1 infectee

  # Assign infectees for Case 1
  infectees <- seq(next_available, by = 1, length.out = R_1)
  infectees <- infectees[infectees > 1 & infectees <= n_cases]  # Ensure valid assignments

  if (length(infectees) > 0) {
    tree$source[infectees] <- tree$id[1]
    tree$R[1] <- length(infectees)
    next_available <- max(infectees) + 1
  }

  # Iterate over remaining cases
  for (i in 2:n_cases) {
    if (next_available > n_cases) break  # Stop if no more cases available

    # Determine how many cases this one infects
    max_possible <- n_cases - next_available + 1
    n_secondary <- min(rpois(1, R_value), max_possible)

    # Identify infectees, ensuring they are AFTER the current case
    infectees <- seq(next_available, by = 1, length.out = n_secondary)
    infectees <- infectees[infectees > i & infectees <= n_cases]  # Ensure no self-infection

    if (length(infectees) > 0) {
      tree$source[infectees] <- tree$id[i]
      tree$R[i] <- length(infectees)
      next_available <- max(infectees) + 1  # Move next available case forward
    }
  }

  # **Ensure the graph is fully connected**
  # Identify cases that were never assigned an infector
  unassigned <- which(is.na(tree$source))

  if (length(unassigned) > 0) {
    for (j in unassigned) {
      # Assign them a random existing case as source
      valid_sources <- tree$id[as.numeric(tree$id) < as.numeric(tree$id[j])]  # Only earlier cases
      if (length(valid_sources) > 0) {
        tree$source[j] <- sample(valid_sources, 1)
      }
    }
  }

  tree$date_onset <- seq_len(n_cases)
  return(tree)
}


test_branching_logic <- function(tree) {
  # Convert id and source to numeric for ordering check
  tree$id <- as.numeric(tree$id)
  tree$source <- as.numeric(tree$source)

  # Remove root cases (NA sources)
  valid_cases <- !is.na(tree$source)

  # Check 1: Every case must be infected by an earlier case
  valid_order <- tree$id[valid_cases] > tree$source[valid_cases]

  # Check 2: Every source ID must exist in the dataset before being referenced
  valid_sources <- tree$source[valid_cases] %in% tree$id[seq_len(sum(valid_cases))]

  # Report test results
  if (all(valid_order) && all(valid_sources)) {
    message("✅ The tree respects the branching logic.")
    return(TRUE)
  } else {
    message("❌ The tree violates branching logic!")
    if (!all(valid_order)) message("   - Some cases are infected by future cases.")
    if (!all(valid_sources)) message("   - Some sources are missing or incorrectly referenced.")
    return(FALSE)
  }
}

# bad_tree_1 <- data.frame(
#   id = as.character(1:10),
#   source = c(NA, 1, 6, 3, 4, 5, 3, 7, 8, 9),  # ❌ Case 3 → Case 6, but Case 6 → Case 3 (cycle)
#   R = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0),
#   date_onset = 1:10
# )
#
# bad_tree_1 <- data.frame(
#   id = as.character(1:10),
#   source = c(NA, 1, 2, 3, 4, 10, 6, 7, 8, 9),  # ❌ Case 5 is infected by Case 10
#   R = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0),
#   date_onset = 1:10
# )
#
# test_branching_logic(bad_tree_1)
#
# bad_tree_1 <- epicontacts::make_epicontacts(
#   linelist = subset(bad_tree_1, select = -source),
#   contacts = subset(bad_tree_1, select = c(source, id)),
#   na_rm_contacts = TRUE,
#   directed = TRUE
# )


####################################
# make TT
####################################
make_tt <- function(n_cases, superspreaders, lambda) {
  if (!all(superspreaders %in% seq_len(n_cases))) {
    stop("Superspreader IDs must be within the range of case IDs.")
  }

  # Initialise the tree with a linear chain
  tree <- data.frame(
    id = as.character(seq_len(n_cases)),
    source = rep(NA, n_cases),
    R = c(rep(1, n_cases-1), 0)
  )

  # Assign sources based on the initial linear chain
  for (i in 2:n_cases) {
    tree$source[i] <- as.character(i - 1)
  }

  # Handle superspreaders
  for (ss in superspreaders) {
    # Determine how many additional cases this superspreader will infect
    additional_cases <- rpois(1, lambda)

    # Find available cases to be reassigned
    possible_targets <- which(tree$source %in% tree$id[ss:n_cases])

    if (length(possible_targets) > 0) {
      reassign_count <- min(additional_cases, length(possible_targets))

      # Reassign sources to the superspreader in their initial order
      for (j in seq_len(reassign_count)) {
        tree$source[possible_targets[j]] <- as.character(ss)
      }
    }
  }
  tree$date_onset <- 1:n_cases
  return(tree)
}



make_linear_tree <- function(n_cases) {
  tree <- data.frame(
    id = as.character(seq_len(n_cases)),
    source = c(NA, as.character(seq_len(n_cases - 1))),
    R = c(rep(1, n_cases - 1), 0),
    date_onset = seq_len(n_cases)
  )

  return(tree)
}




#
#
# # Test by tree topologies
#
# Gradual Introduction of super-spreading events to control for changing tree topologies.
#
# We start with a linear tree (`A`) and introduce superspreading events or increased branching sizes to create epidemiologically sensible variations in tree topologies.
#
# ## Topology A
#
# A linear chain of infected individuals.
#
# ```{r}
# source(here::here("analysis", "performance", "R", "make_tree.R"))
# A <- make_linear_tree(n_cases = 20)
# A %>% make_epicontacts %>% plot_epicontacts
# ```
#
# ## Topology B
#
# ```{r}
# B <- make_tree(n_cases = 20, R_value = 2)
# B %>% make_epicontacts %>% plot_epicontacts
# ```
#
# ## Topology C
#
# ```{r}
# C <- make_tree(n_cases = 20, R_value = 10)
# C %>% make_epicontacts %>% plot_epicontacts
# ```
#
# ## Topology D
#
# ```{r}
#
# D <- make_tree(n_cases = 20, R_value = 100)
# D %>% make_epicontacts %>% plot_epicontacts
# ```
#
# ## Alternatively add individual super-spreaders
#
# ```{r}
# #set.seed(123)
# n_cases <- 20
# superspreaders <- c(18)
# lambda <- 100  #R for superspreaders
# make_tt(n_cases, superspreaders, lambda) %>% make_epicontacts() %>%  plot_epicontacts()
# ```
#
# ```{r}
# set.seed(123)
# n_cases <- 20
# superspreaders <- c(17)
# lambda <- 100
# make_tt(n_cases, superspreaders, lambda) %>% make_epicontacts() %>%  plot_epicontacts()
# ```
#
# ```{r}
# set.seed(123)
# n_cases <- 20
# superspreaders <- c(3, 9)
# lambda <- 5
# make_tt(n_cases, superspreaders, lambda) %>% make_epicontacts() %>%  plot_epicontacts()
# ```
#
