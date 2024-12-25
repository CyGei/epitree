# Euclidean distance between two patristic matrices
euclidean <- function(mat1, mat2) {
  sqrt(sum((mat1[lower.tri(mat1)] - mat2[lower.tri(mat2)]) ^ 2))
}
