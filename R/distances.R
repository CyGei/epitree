# Compute the patristic distance matrix for a single tree
patristic <- function(tree) {
  g <- graph_from_data_frame(tree, directed = TRUE)
  igraph::distances(g, mode = "all")
}

# Euclidean distance between two patristic matrices
euclidean <- function(mat1, mat2) {
  sqrt(sum((mat1 - mat2) ^ 2))
}


# Compare chains of posterior trees
compare_chains <- function(...,
                           f_within = patristic,
                           f_between = euclidean,
                           adonis2_args = list()) {
  chains <- list(...)

  # Check for valid input
  if (length(chains) < 2)
    stop("At least two chains are required.")
  if (is.null(names(chains)))
    names(chains) <- as.character(seq_along(chains))

  # Generate chain_id for each data frame in the chains
  chain_lengths <- sapply(chains, function(chain)
    ifelse(is.data.frame(chain), 1, length(chain)))
  chain_id <- factor(rep(names(chains), chain_lengths))

  if (all(sapply(chains, is.data.frame))) {
    master_chain <- chains # if only 1 dataframe per chain
  } else {
    master_chain <- unlist(chains, recursive = FALSE)
  }

  within_list <- lapply(master_chain, f_within)

  # Get all unique pair combinations#
  n <- length(master_chain)
  pairs <- t(combn(n, 2))
  d <- matrix(0, n, n)


  distances <- apply(pairs, 1, function(p) {
    f_between(within_list[[p[1]]], within_list[[p[2]]])
  })
  d[lower.tri(d)] <- distances

  # adonis2
  if (!is.list(adonis2_args)) {
    stop("adonis2_args must be a list.")
  }
  if (length(adonis2_args) < 1) {
    adonis <- vegan::adonis2(formula = stats::as.dist(d) ~ chain_id,
                             data = data.frame(chain_id = chain_id))
  } else {
    # TODO: Ask @Thibaut why I get that warning message, I suspect its the formula not passing well
    adonis <- do.call(vegan::adonis2, c(
      list(
        formula = stats::as.dist(d) ~ chain_id,
        data = data.frame(chain_id = chain_id)
      ),
      adonis2_args
    ))
  }
  return(adonis)
}




get_chisq <- function(chain_x, chain_y, levels = NULL) {
  make_ttab <- function(chain, levels) {
    tab <- do.call(rbind, chain)
    tab <- o2ools::ttable(from = tab$from,
                          to = tab$to,
                          levels = levels)
    tab <- as.data.frame(tab)
    tab <- tab[tab$from != tab$to, ]
    return(tab)
  }

  tab_x <- make_ttab(chain_x, levels)
  tab_y <- make_ttab(chain_y, levels)

  tab <- merge(tab_x, tab_y, by = c("from", "to"), all = TRUE)
  tab <- tab[!(tab$Freq.x == 0 & tab$Freq.y == 0), ]

  return(chisq.test(tab$Freq.x, tab$Freq.y))
}
