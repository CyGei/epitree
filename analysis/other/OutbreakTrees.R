# https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001685#pbio.3001685
devtools::load_all()
pacman::p_load(tidyverse, httr, igraph)
# URL of the RDS file on GitHub
url <- "https://github.com/DrakeLab/taube-transmission-trees/raw/master/data/data_tibble_trees.RDS"
data <- readRDS(url(url, "rb"))
trees <- data %>%
  filter(Size >= 30) %>%
  .[["tree"]] %>%
  lapply(., \(g) igraph::upgrade_graph(g) |> igraph::as_long_data_frame() |> select(from, to))
invisible(lapply(trees, check_tree))
sapply(trees, nrow) |> DescTools::Mode()

N <- length(trees)
n_each <- N %/% 2
idx <- sample(1:N, n_each, replace = FALSE) #idx2 <- setdiff(1:N, idx1)
chain1 <- trees[idx]
chain2 <- trees[-idx]
compare_chains(chain1, chain2)


