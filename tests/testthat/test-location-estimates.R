test_that("spatial.median minimizes the sum of distances and is orthogonally equivariant", {
  set.seed(21)
  X <- rmvnorm(50, c(1, 2, 3), diag(c(1, 2, 3)))
  m <- spatial.median(X, eps = 1e-10)
  obj <- function(mu) sum(sqrt(colSums((t(X) - mu)^2)))
  for (k in 1:10) expect_lte(obj(m), obj(m + rnorm(3, sd = 0.01)))

  Q <- qr.Q(qr(matrix(rnorm(9), 3)))
  mQ <- spatial.median(X %*% t(Q) + matrix(b3, 50, 3, byrow = TRUE), eps = 1e-10)
  expect_equal(mQ, drop(Q %*% m + b3), tolerance = 1e-6)

  # univariate data give the median
  expect_equal(spatial.median(matrix(c(3, 1, 2, 10, 7))), 3)
})

test_that("spatial.median of centrally symmetric data is the center", {
  set.seed(22)
  Y <- rmvnorm(20, c(0, 0, 0))
  X <- rbind(Y, -Y) + matrix(b3, 40, 3, byrow = TRUE)
  expect_equal(spatial.median(X, eps = 1e-10), b3, tolerance = 1e-6)
})

test_that("HR.Mest is affine equivariant and solves its estimating equations", {
  set.seed(23)
  X <- rmvnorm(60, c(0, 0, 0), diag(c(1, 2, 3)))
  hr <- HR.Mest(X, eps.scale = 1e-8, eps.center = 1e-8)
  hrA <- HR.Mest(X %*% t(A3) + matrix(b3, 60, 3, byrow = TRUE),
                 eps.scale = 1e-8, eps.center = 1e-8)
  expect_equal(unname(hrA$center), drop(A3 %*% hr$center + b3), tolerance = 1e-5)
  expect_equal(shp(hrA$scatter), shp(A3 %*% hr$scatter %*% t(A3)), tolerance = 1e-5)

  U <- spatial.sign(X, center = hr$center, shape = hr$scatter)
  expect_equal(colMeans(U), rep(0, 3), tolerance = 1e-4)
  expect_equal(3 * crossprod(U) / 60, diag(3), tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("hl.loc equals the median of the Walsh averages", {
  set.seed(24)
  x <- rnorm(15)
  expect_equal(hl.loc(x), unname(wilcox.test(x, conf.int = TRUE)$estimate))
})

# brute-force reference for vdw.loc: evaluate the van der Waerden signed-rank
# statistic on every interval between neighbouring Walsh averages and take
# (sup{T > 0} + inf{T < 0}) / 2
vdw.ref <- function(x) {
  n <- length(x)
  s <- outer(x, x, "+") / 2
  w <- sort(unique(s[upper.tri(s, diag = TRUE)]))
  if (length(w) == 1) return(w)
  mids <- (w[-1] + w[-length(w)]) / 2
  Tm <- sapply(mids, function(mu) {
    z <- x - mu
    sum(sign(z) * qnorm(0.5 + rank(abs(z)) / (2 * (n + 1))))
  })
  (w[sum(Tm > 1e-9) + 1] + w[sum(Tm >= -1e-9) + 1]) / 2
}

test_that("vdw.loc agrees with a brute-force evaluation of its definition", {
  set.seed(25)
  for (i in 1:100) {
    x <- rt(sample(1:40, 1), 3)
    expect_equal(vdw.loc(x), vdw.ref(x))
  }
  # data with ties
  for (i in 1:50) {
    x <- round(rnorm(sample(2:30, 1)), 1)
    expect_equal(vdw.loc(x), vdw.ref(x))
  }
})

test_that("vdw.loc is exactly equivariant and recovers the center of symmetric samples", {
  set.seed(26)
  x <- rt(30, 4)
  expect_equal(vdw.loc(3 * x + 5), 3 * vdw.loc(x) + 5)
  expect_equal(vdw.loc(-x), -vdw.loc(x))
  y <- abs(rnorm(15))
  expect_equal(vdw.loc(c(2 - y, 2 + y)), 2)
  expect_equal(vdw.loc(c(2 - y, 2, 2 + y)), 2)
})

test_that("vdw.loc handles small samples and gives finite values", {
  expect_equal(vdw.loc(3), 3)
  expect_equal(vdw.loc(c(1, 4)), 2.5)
  expect_equal(vdw.loc(c(5, 5, 5)), 5)
  expect_error(vdw.loc(numeric(0)), "at least one")
  # up to version 1.1-3 about 0.4 percent of such samples gave -Inf or Inf
  set.seed(3)
  est <- replicate(100, vdw.loc(rt(sample(5:200, 1), 3)))
  expect_true(all(is.finite(est)))
  # int.diff and maxiter are no longer used
  x <- rnorm(20)
  expect_identical(vdw.loc(x, int.diff = 2, maxiter = 1), vdw.loc(x))
  expect_equal(vdw.loc(c(x, NA), na.action = na.omit), vdw.loc(x))
})
