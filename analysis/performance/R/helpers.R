
# Outbreak Simulation -----------------------------------------------------
#' Simulate outbreaks within a target size range
#'
#' This function repeatedly simulates outbreaks until a specified number of simulations
#' is obtained, ensuring each simulation's outbreak size falls within a tolerance range around
#' a target size.
#'
#' @param target_size Integer. The target outbreak size.
#' @param R_values Integer vector. Reproduction number(s) used in the simulation.
#' @param tolerance Numeric. Tolerance for outbreak size (default is 0.2 for ±20% of target size).
#' @param n_simulations Integer. Number of valid simulations to generate (default is 200).
#' @param max_attempts Integer. Maximum number of simulation attempts before checking progress (default is 1000).
#' @param try_gain Numeric. Fraction of target simulations required to extend max_attempts (default is 1, i.e. no extension).
#'
#' @return A list of outbreak simulations (data frames) whose sizes are within
#'   \code{target_size * (1 - tolerance)} and \code{target_size * (1 + tolerance)}.
#'
#' @details The function simulates outbreaks using \code{simulacr::simulate_outbreak} until
#'   \code{n_simulations} valid simulations are collected. If the number of attempts reaches
#'   \code{max_attempts} and at least \code{try_gain * n_simulations} simulations have been obtained,
#'   an additional 1000 attempts are added automatically; otherwise, the function issues a warning.
#'
#' @examples
#' \dontrun{
#'   # Simulate 100 outbreaks targeting 50 cases ±20% using negative binomial R_values.
#'   sims <- simulate_outbreaks(target_size = 50,
#'                              R_values = rnbinom(100, size = 0.2, mu = 3),
#'                              tolerance = 0.2,
#'                              n_simulations = 100,
#'                              max_attempts = 1000,
#'                              try_gain = 0.8)
#' }
#'
simulate_outbreaks <- function(target_size,
                               R_values,
                               tolerance = 0.2,
                               n_simulations = 200,
                               max_attempts = 1000,
                               try_gain = 1) {
  # Calculate acceptable size range
  min_size <- round(target_size * (1 - tolerance))
  max_size <- round(target_size * (1 + tolerance))
  sims <- list()
  attempts <- 0

  while (length(sims) < n_simulations) {
    # Check if current max_attempts limit is reached
    if (attempts >= max_attempts) {
      if (length(sims) >= try_gain * n_simulations) {
        message(
          "Reached ",
          try_gain * 100,
          "% of target simulations. Adding an additional 1000 attempts."
        )
        max_attempts <- max_attempts + 1000
      } else {
        warning(
          "Maximum number of attempts reached at ",
          length(sims),
          "simulations."
        )
        break
      }
    }

    sim <- simulacr::simulate_outbreak(
      duration = 100,
      population_size = round(target_size / epitrix::R02AR(mean(R_values))),
      R_values = R_values,
      dist_incubation = outbreaker2::fake_outbreak$w,
      dist_generation_time = outbreaker2::fake_outbreak$w
    )$data %>%
      relabel_tree(id_cols = c("id", "source"), date_col = "date_onset")
    sim <- shift_init_date(sim)

    attempts <- attempts + 1

    # Accept simulation if within acceptable size range
    if (nrow(sim) >= min_size && nrow(sim) <= max_size) {
      sims[[length(sims) + 1]] <- sim
    }
  }

  return(sims)
}


# Process Simulations -----------------------------------------------------
#' @title Match data frames by row counts
#'
#' @description Matches and pairs data frames from lists \code{A} and \code{B} based on identical row counts.
#'
#' @param A A list of data frames.
#' @param B A list of data frames.
#' @param n Number of pairs to return.
#'
#' @return A list of length \code{n}, each element containing a pair of data frames with the same number of rows.
matching_pairs <- function(A, B, n = 100) {
  nrows_A <- sapply(A, nrow)
  nrows_B <- sapply(B, nrow)

  dfA <- data.frame(index_A = seq_along(A), nrows = nrows_A)
  dfB <- data.frame(index_B = seq_along(B), nrows = nrows_B)

  merged_df <- merge(dfA, dfB, by = "nrows")

  #keep unique index from A @ASK if necessary
  #merged_df <- merged_df[!duplicated(merged_df$index_A), ]


  if (nrow(merged_df) > n) {
    merged_df <- merged_df[sample(1:nrow(merged_df), n), ]
  }

  paired_list <- mapply(function(a_idx, b_idx) {
    list(A = A[[a_idx]], B = B[[b_idx]])
  }, merged_df$index_A, merged_df$index_B, SIMPLIFY = FALSE)

  if (length(paired_list) < n) {
    stop("Not enough matching pairs found to create ", n, "pairs.")
  }

  # Check that each pair has the same number of rows
  rows_match <- sapply(paired_list, function(pair) nrow(pair$A) == nrow(pair$B))
  if (!all(rows_match)) {
    stop("Not all pairs have matching row counts.")
  }

  return(paired_list)
}


# Add a day for cases with the same onset date as the first case
# This avoids having multiple introductions for outbreaker2
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


# Outbreak Reconstruction -------------------------------------------------
# Helper function to run outbreaker on a single simulation
run_outbreaker <- function(sim, ctd_fraction = 0.5) {
  n_cases <- nrow(sim)
  # Adding ctd data to improve reconstruction
  ctd_sample <- sample(2:n_cases, round(n_cases * ctd_fraction))

  config <- outbreaker2::create_config(
    init_tree = "star",
    move_mu = FALSE,
    init_pi = 1,
    move_pi = FALSE,
    init_kappa = 1,
    move_kappa = FALSE,
    find_import = FALSE
  )

  data <- outbreaker2::outbreaker_data(
    dates = sim$date_onset,
    w_dens = outbreaker2::fake_outbreak$w,
    f_dens = outbreaker2::fake_outbreak$w,
    ctd = sim[ctd_sample, c("source", "id")],
    ids = sim$id
  )

  out <- outbreaker2::outbreaker(data = data, config = config)
  out <- o2ools::identify(out, sim$id)
  return(out)
}
