#' Calculate the Euclidean distance between two distance matrices.
#'
#' This function computes the Euclidean distance between the lower triangular 
#' parts of two given matrices. 
#'
#' @param mat1 A numeric matrix.
#' @param mat2 A numeric matrix.
#'
#' @return A numeric value representing the Euclidean distance between the 
#' lower triangular parts of \code{mat1} and \code{mat2}.
#'
#' @examples
#' mat1 <- matrix(c(1, 2, 3, 4), 2, 2)
#' mat2 <- matrix(c(4, 3, 2, 1), 2, 2)
#' euclidean(mat1, mat2)
#' @export
euclidean <- function(mat1, mat2) {
  sqrt(sum((mat1[lower.tri(mat1)] - mat2[lower.tri(mat2)]) ^ 2))
}
