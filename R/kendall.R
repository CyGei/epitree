pacman::p_load(treespace)
## a simple who infected whom matrix:
tree1 <- cbind(from = 1:5, to = 2:6)
tree2 <- cbind(from = c(1, 5, 2, 2, 3), to = 2:6)
tree3 <- cbind(from = c(2, 2, 2, 2, 6), to = c(1, 3, 4, 6, 5))
m1 <- findMRCIs(tree1) # find the source case, MRCIs and MRCI depths for tree 1
m2 <- findMRCIs(tree2)
m3 <- findMRCIs(tree3)
matList <- list(m1$mrciDepths, m2$mrciDepths, m3$mrciDepths) # create a list of the mrciDepths matrices
wiwTreeDist(matList, sampled = 1:6)
euclidean(matList[[1]], matList[[2]])
pacman::p_load(outbreaker2, ape, igraph, vegan, distcrete, ggplot2, scales)

chains = list(
  as.data.frame(tree1),
  as.data.frame(tree2),
  as.data.frame(tree3)
)

chains = list(
  list(
  as.data.frame(tree1),
  as.data.frame(tree2),
  as.data.frame(tree3)
),
  list(
  as.data.frame(tree1),
  as.data.frame(tree2),
  as.data.frame(tree3)
)
)

compare_chains()
mrciDepth <- function(tree) {
  treespace::findMRCIs(as.matrix(tree))$mrciDepths
}

# > vegan::adonis2(formula = stats::as.dist(d) ~ chain_id,
#                  +                              data = data.frame(chain_id = chain_id))
# 'nperm' >= set of all permutations: complete enumeration.
# Set of permutations < 'minperm'. Generating entire set.
# No residual component
#
# vegan::adonis2(formula = stats::as.dist(d) ~ chain_id, data = data.frame(chain_id = chain_id))
# Df SumOfSqs R2 F Pr(>F)
# Model     2   68.667  1
# Residual  0    0.000  0
# Total     2   68.667  1
