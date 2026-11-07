# Internal helpers for the compiled code and the spatial sign functions

# coerce to a numeric matrix with double storage, as expected by src/icsnp.c
num.matrix <- function(X)
    {
    X <- as.matrix(X)
    if (!is.numeric(X)) stop("'X' must be numeric")
    storage.mode(X) <- "double"
    X
    }

# Euclidean norms of the rows of X
row.norms <- function(X) sqrt(rowSums(X^2))

# sum of the outer products of the spatial signs of the rows of X
sumsignout <- function(X) crossprod(X / row.norms(X))
