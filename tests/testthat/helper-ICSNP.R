# Helpers used in several test files

# standardize a scatter matrix to have determinant one
shp <- function(S) S / det(S)^(1 / ncol(S))

# a fixed nonsingular, non-orthogonal transformation matrix and shift
A3 <- matrix(c(2, 1, 0.5, -1, 3, 0, 0.5, 0.2, 1.5), ncol = 3)
b3 <- c(1, -2, 3)

# reference implementation of all pairwise differences, sums, products
# (rows ordered as (1,2), (1,3), ..., (1,n), (2,3), ...)
pair.ref <- function(X, FUN) {
  idx <- utils::combn(nrow(X), 2)
  FUN(X[idx[1, ], , drop = FALSE], X[idx[2, ], , drop = FALSE])
}
