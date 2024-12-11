# Compute the patristic distance matrix for a single tree
patristic <- function(tree) {
  g <- graph_from_data_frame(tree, directed = TRUE)
  igraph::distances(g, mode = "all")
}

# Euclidean distance between two patristic matrices
euclid <- function(mat1, mat2) {
  sqrt(sum((mat1 - mat2)^2))
}


# Compare chains of posterior trees
compare_chains <- function(..., f = euclid, use_for = FALSE, adonis2_args = list()) {
  # browser()
  chains <- list(...)

  # Check for valid input
  if (length(chains) < 2) stop("At least two chains are required.")
  if (is.null(names(chains))) names(chains) <- as.character(seq_along(chains))

  chain_id <- factor(rep(names(chains), vapply(chains, length, integer(1))))

  master_chain <- unlist(chains, recursive = FALSE)

  patristic_list <- lapply(master_chain, patristic)

  # Get all unique pair combinations#
  n <- length(master_chain)
  pairs <- t(combn(n, 2))
  d <- matrix(0, n, n)

  if (use_for == TRUE) {
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        d[j, i] <- f(patristic_list[[i]], patristic_list[[j]])
      }
    }
  } else {
    distances <- apply(pairs, 1, function(p) {
      f(patristic_list[[p[1]]], patristic_list[[p[2]]])
    })
    d[lower.tri(d)] <- distances
  }

  # adonis2
  if (!is.list(adonis2_args)) {
    stop("adonis2_args must be a list.")
  }
  if (length(adonis2_args) < 1) {
    adonis <- vegan::adonis2(
      formula = stats::as.dist(d) ~ chain_id,
      data = data.frame(chain_id = chain_id)
    )
  } else { #TODO: Ask @Thibaut why I get that warning message, I suspect its the formula not passing well
    adonis <- do.call(
      vegan::adonis2,
      c(
        list(
          formula = stats::as.dist(d) ~ chain_id,
          data = data.frame(chain_id = chain_id)
        ),
        adonis2_args
      )
    )
  }
  return(adonis)
}




get_chisq <- function(chain_x, chain_y, levels = NULL) {
  make_ttab <- function(chain, levels) {
    tab <- do.call(rbind, chain)
    tab <- o2ools::ttable(from = tab$from, to = tab$to, levels = levels)
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



# #######################################################
# # Initial Implementation SLOW!!
# #######################################################
# # Compare trees ----------------------------------------------------

# compare_trees <- function(tree1, tree2) {
#   euclid <- function(mat1, mat2) {
#     sqrt(sum((mat1 - mat2)^2))
#   }

#   # convert to graphs to compute patristic ('generational') distance
#   g1 <- igraph::graph_from_data_frame(tree1, directed = TRUE)
#   g2 <- igraph::graph_from_data_frame(tree2, directed = TRUE)

#   # computes the euclidian distance between 2 patristic matrices
#   euclid(
#     igraph::distances(g1, mode = "all"),
#     igraph::distances(g2, mode = "all")
#   )
# }



# # Compare chains of posterior trees ----------------------------------------------------

# compare_chains <- function(..., f = compare_trees, ...) {
#   chains <- list(...)

#   if (length(chains) < 2) {
#     stop("At least two chains are required for comparison.")
#   }

#   master_chain <- unlist(chains, recursive = FALSE)

#   chain_id <- data.frame(chaind_id = rep(names(chains), sapply(chains, length)))

#   d <- as.dist(outer(master_chain, master_chain, Vectorize(f)))

#   result <- vegan::adonis2(d ~ chain_id, data = chain_id, ...)

#   return(result)
# }

# # compare_chains <- function(chain1, chain2, ...) {
# #   master_chain <- c(chain1, chain2)
# #   chain_id <- data.frame(
# #     chain_id = rep(c("chain1", "chain2"),
# #       times = c(length(chain1), length(chain2))
# #     )
# #   )
# #   d <- as.dist(outer(master_chain, master_chain, Vectorize(compare_trees)))
# #   vegan::adonis2(d ~ chain_id, data = chain_id, ...)
# # }

# set.seed(123)
# chain1 = sample(trees, size = 100, replace = TRUE)
# chain2 = sample(trees, size = 10, replace = TRUE)
# compare_chains(chain1, chain2)
