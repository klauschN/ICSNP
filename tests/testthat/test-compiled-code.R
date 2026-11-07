# These tests compare each routine in src/ with a plain R reference
# implementation, so that the compiled code can be refactored safely.

test_that("pair.diff, pair.sum and pair.prod agree with R references", {
  set.seed(41)
  X <- matrix(rnorm(7 * 3), 7)
  expect_equal(pair.diff(X), pair.ref(X, `-`), ignore_attr = TRUE)
  expect_equal(pair.sum(X), pair.ref(X, `+`), ignore_attr = TRUE)
  expect_equal(pair.prod(X), pair.ref(X, `*`), ignore_attr = TRUE)
  expect_equal(dim(pair.diff(X)), c(choose(7, 2), 3))
  # a single column
  x <- matrix(rnorm(5))
  expect_equal(pair.diff(x), pair.ref(x, `-`), ignore_attr = TRUE)
  # integer input, data frames and fewer than two rows
  Xi <- matrix(1:12, 4)
  expect_equal(pair.prod(Xi), pair.ref(Xi, `*`), ignore_attr = TRUE)
  expect_equal(pair.sum(as.data.frame(X)), pair.sum(X))
  expect_equal(dim(pair.diff(X[1, , drop = FALSE])), c(0, 3))
  expect_error(pair.diff(matrix(letters[1:4], 2)), "numeric")
  # more pairs than a matrix can have rows gives an error before allocation
  expect_error(pair.diff(matrix(0, 65537, 1)), "too many pairs")
})

test_that("row norms and sums of sign outer products agree with R references", {
  set.seed(42)
  X <- matrix(rnorm(20 * 3), 20)
  expect_equal(ICSNP:::row.norms(X), sqrt(rowSums(X^2)))
  U <- X / sqrt(rowSums(X^2))
  expect_equal(ICSNP:::sumsignout(X), crossprod(U), ignore_attr = TRUE)
})

test_that("sum of pairwise difference sign outer products agrees with R reference", {
  set.seed(43)
  X <- matrix(rnorm(15 * 3), 15)
  D <- pair.ref(X, `-`)
  U <- D / sqrt(rowSums(D^2))
  expect_equal(ICSNP:::SSCov(X), crossprod(U) / nrow(D), ignore_attr = TRUE)
  expect_equal(attr(ICSNP:::SSCov(X), "n.zero"), 0)
  # tied rows: the zero differences are skipped and counted
  Xt <- rbind(X, X[1, ], X[2, ], X[2, ])
  Dt <- pair.ref(Xt, `-`)
  Dt <- Dt[rowSums(Dt^2) > 0, ]
  Ut <- Dt / sqrt(rowSums(Dt^2))
  res <- ICSNP:::SSCov(Xt)
  expect_equal(attr(res, "n.zero"), 4)
  expect_equal(res, crossprod(Ut) / choose(nrow(Xt), 2), ignore_attr = TRUE)
})

test_that("one symmetrized Huber step agrees with R reference", {
  set.seed(44)
  X <- matrix(rnorm(15 * 3), 15)
  V0 <- cov(X)
  p <- 3; qg <- 0.9
  cs <- 2 * qchisq(qg, p)
  sig <- 2 * pchisq(cs / 2, p + 2) + (cs / p) * (1 - qg)
  D <- pair.ref(X, `-`)
  d2 <- mahalanobis(D, 0, V0)
  w <- ifelse(d2 <= cs, 1 / sig, (cs / d2) / sig)
  ref <- crossprod(sqrt(w) * D) / nrow(D)
  expect_equal(ICSNP:::SSCov.hub(X, V0, cs, sig), ref, ignore_attr = TRUE)
})

# pure R versions of the iterations, as done in the R code up to 1.1-3
# (with symmetric square roots instead of the Cholesky factors used in C)
duembgen.ref <- function(X, V, n.steps) {
  D <- pair.ref(X, `-`)
  D <- D[rowSums(D^2) > 0, , drop = FALSE]
  n.pairs <- choose(nrow(X), 2)
  for (s in seq_len(n.steps)) {
    e <- eigen(V, symmetric = TRUE)
    sq <- e$vectors %*% (sqrt(e$values) * t(e$vectors))
    E <- D %*% solve(sq)
    S <- crossprod(E / sqrt(rowSums(E^2))) / n.pairs
    V <- shp(sq %*% S %*% sq)
  }
  V
}
huber.ref <- function(X, V, n.steps, qg = 0.9) {
  p <- ncol(X)
  cs <- 2 * qchisq(qg, p)
  sig <- 2 * pchisq(cs / 2, p + 2) + (cs / p) * (1 - qg)
  D <- pair.ref(X, `-`)
  for (s in seq_len(n.steps)) {
    d2 <- mahalanobis(D, 0, V)
    w <- ifelse(d2 <= cs, 1 / sig, (cs / d2) / sig)
    V <- crossprod(sqrt(w) * D) / nrow(D)
  }
  V
}

test_that("the iterations in C agree with R references step by step", {
  set.seed(45)
  X <- rmvnorm(25, c(0, 0, 0), diag(c(1, 4, 9)))
  V0 <- shp(cov(X))
  for (s in c(1, 2, 5)) {
    expect_equal(duembgen.shape(X, steps = s), duembgen.ref(X, V0, s),
                 tolerance = 1e-10, ignore_attr = TRUE)
  }
  # with ties
  Xt <- rbind(X, X[1:2, ])
  expect_equal(suppressWarnings(duembgen.shape(Xt, steps = 3)),
               duembgen.ref(Xt, shp(cov(Xt)), 3), tolerance = 1e-10, ignore_attr = TRUE)
  # the converged symmetrized Huber estimate is a fixed point of the R step,
  # and agrees with a fixed number of R steps
  sh <- symm.huber(X, eps = 1e-12)
  expect_equal(huber.ref(X, sh, 1), sh, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(huber.ref(X, cov(X), 500), sh, tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("the iterations in C report maxiter and non positive definite matrices", {
  set.seed(46)
  X <- rmvnorm(30, c(0, 0))
  expect_error(duembgen.shape(X, maxiter = 1, eps = 1e-12), "maxiter reached")
  expect_error(symm.huber(X, maxiter = 1, eps = 1e-12), "maxiter reached without convergence")
  expect_error(symm.huber(X, init = diag(c(1, -1))), "not positive definite")
  expect_error(duembgen.shape(X, init = diag(c(1, 0))), "positive definite|infinite|NaN")
  # the result is a plain matrix with the column names of X
  colnames(X) <- c("a", "b")
  expect_identical(names(attributes(duembgen.shape(X))), c("dim", "dimnames"))
  expect_identical(colnames(symm.huber(X)), c("a", "b"))
})

test_that("spatial.sign returns unit vectors and zero for points at the center", {
  X <- rbind(c(0, 0), c(1, 1), c(2, 0), c(-3, 4))
  U <- spatial.sign(X, center = c(0, 0), shape = FALSE)
  expect_equal(unname(U[1, ]), c(0, 0))
  expect_equal(sqrt(rowSums(U[-1, ]^2)), rep(1, 3), ignore_attr = TRUE)
  expect_equal(attr(U, "center"), c(0, 0))
  expect_equal(drop(spatial.sign(matrix(c(-2, 0, 3)), center = 0)), c(-1, 0, 1),
               ignore_attr = TRUE)
})
