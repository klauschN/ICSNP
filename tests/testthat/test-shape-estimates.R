test_that("tyler.shape has determinant one and solves its estimating equation", {
  set.seed(31)
  X <- rmvnorm(80, c(0, 0, 0), diag(c(1, 2, 3)))
  V <- tyler.shape(X, location = c(0, 0, 0), eps = 1e-10)
  expect_equal(det(V), 1, tolerance = 1e-8)
  expect_equal(V, t(V))
  U <- spatial.sign(X, center = c(0, 0, 0), shape = V)
  expect_equal(3 * crossprod(U) / 80, diag(3), tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("tyler.shape is affine equivariant and the C and R versions agree", {
  set.seed(32)
  X <- rmvnorm(80, c(0, 0, 0), diag(c(1, 2, 3)))
  V <- tyler.shape(X, location = c(0, 0, 0), eps = 1e-10)
  VA <- tyler.shape(X %*% t(A3) + matrix(b3, 80, 3, byrow = TRUE), location = b3, eps = 1e-10)
  expect_equal(VA, shp(A3 %*% V %*% t(A3)), tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(tyler.shape(X, location = c(0, 0, 0), eps = 1e-10, in.R = TRUE), V,
               tolerance = 1e-6)
})

test_that("duembgen.shape: C and R versions agree, affine equivariance, no location needed", {
  set.seed(33)
  X <- rmvnorm(60, c(0, 0, 0), diag(c(1, 2, 3)))
  D <- duembgen.shape(X, eps = 1e-10)
  expect_equal(det(D), 1, tolerance = 1e-8)
  expect_equal(duembgen.shape(X, eps = 1e-10, in.R = TRUE), D, tolerance = 1e-6)
  DA <- duembgen.shape(X %*% t(A3) + matrix(b3, 60, 3, byrow = TRUE), eps = 1e-10)
  expect_equal(DA, shp(A3 %*% D %*% t(A3)), tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("duembgen.shape handles tied observations (bug fixed in 1.2-0)", {
  set.seed(34)
  X <- rmvnorm(30, c(0, 0, 0))
  X <- rbind(X, X[1, ], X[2, ])
  expect_warning(D <- duembgen.shape(X, eps = 1e-10), "2 pairwise")
  expect_true(all(is.finite(D)))
  expect_warning(D.R <- duembgen.shape(X, eps = 1e-10, in.R = TRUE))
  expect_equal(D, D.R, tolerance = 1e-6)
})

test_that("weighted Duembgen shape matrix reduces to the unweighted one", {
  set.seed(35)
  X <- rmvnorm(40, c(0, 0, 0), diag(c(1, 2, 3)))
  D <- duembgen.shape(X, eps = 1e-10)
  expect_equal(duembgen.shape.wt(X, eps = 1e-10), D, tolerance = 1e-5)
  # zero weights remove observations
  w <- rep(c(1, 0), c(30, 10))
  expect_equal(duembgen.shape.wt(X, w, eps = 1e-10),
               duembgen.shape(X[1:30, ], eps = 1e-10), tolerance = 1e-5)
})

test_that("symm.huber is affine equivariant and agrees with its weighted version", {
  set.seed(36)
  X <- rmvnorm(50, c(0, 0, 0), diag(c(1, 2, 3)))
  S <- symm.huber(X, eps = 1e-10)
  SA <- symm.huber(X %*% t(A3) + matrix(b3, 50, 3, byrow = TRUE), eps = 1e-10)
  expect_equal(SA, A3 %*% S %*% t(A3), tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(symm.huber.wt(X, eps = 1e-10), S, tolerance = 1e-6)
})

test_that("symm.huber is consistent for the covariance matrix at the normal model", {
  set.seed(37)
  Sigma <- matrix(c(2, 0.5, 0.5, 1), 2)
  X <- rmvnorm(2000, c(0, 0), Sigma)
  expect_equal(symm.huber(X), Sigma, tolerance = 0.1, ignore_attr = TRUE)
})

test_that("HP1.shape returns a positive definite shape matrix with determinant one", {
  set.seed(38)
  X <- rmvnorm(50, c(0, 0, 0), diag(c(1, 2, 3)))
  V <- HP1.shape(X)
  expect_equal(det(V), 1, tolerance = 1e-8)
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(all(eigen(V, symmetric = TRUE)$values > 0))
})
