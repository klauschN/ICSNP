test_that("ind.ctest is invariant under increasing marginal transformations", {
  set.seed(11)
  X <- cbind(rnorm(50), rnorm(50), rexp(50), rnorm(50))
  X[, 3] <- X[, 3] + X[, 1]
  Xt <- cbind(exp(X[, 1]), X[, 2]^3, log(X[, 3] + 10), 2 * X[, 4] + 1)
  for (s in c("rank", "normal")) {
    expect_equal(ind.ctest(Xt, 1:2, scores = s)$statistic,
                 ind.ctest(X, 1:2, scores = s)$statistic, tolerance = 1e-10)
  }
  # index2 = NULL means the complement of index1
  expect_equal(ind.ctest(X, 1:2)$statistic, ind.ctest(X, 1:2, 3:4)$statistic)
  expect_equal(unname(ind.ctest(X, 1:2)$parameter), 4)
})

test_that("ind.ictest works when one of the blocks is univariate (bug fixed in 1.2-0)", {
  set.seed(12)
  W <- cbind(rmvnorm(60, rep(0, 3)), rt(60, 4))
  for (s in c("sign", "rank", "normal")) {
    for (m in c("approximation", "permutation", "simulation")) {
      p1 <- ind.ictest(W, 1:3, 4, scores = s, method = m, n.simu = 20)$p.value
      p2 <- ind.ictest(W, 4, 1:3, scores = s, method = m, n.simu = 20)$p.value
      expect_true(p1 >= 0 && p1 <= 1, label = paste(s, m, "second block univariate"))
      expect_true(p2 >= 0 && p2 <= 1, label = paste(s, m, "first block univariate"))
    }
  }
})

test_that("ind.ictest resampled statistics are scaled like the observed one (bug fixed in 1.2-0)", {
  # under independence the statistic is approximately chi-square with p1 * p2
  # degrees of freedom, so the resampled statistics should have mean close to p1 * p2
  set.seed(13)
  n <- 200; p1 <- 2; p2 <- 2
  Z1 <- matrix(rnorm(n * p1), n)
  Z2 <- matrix(rnorm(n * p2), n)
  args <- list(apply(Z1, 2, sign), apply(abs(Z1), 2, rank),
               apply(Z2, 2, sign), apply(abs(Z2), 2, rank), n)
  for (f in c(".Q.perm.sign", ".Q.perm.rank", ".Q.perm.normal")) {
    Q <- replicate(300, do.call(getFromNamespace(f, "ICSNP"), args))
    expect_true(abs(mean(Q) / (p1 * p2) - 1) < 0.25, label = f)
  }
  for (f in c(".Q.simu.sign", ".Q.simu.rank", ".Q.simu.normal")) {
    FUN <- getFromNamespace(f, "ICSNP")
    Q <- replicate(300, FUN(matrix(rnorm(n * (p1 + p2)), n), p1, p2, n))
    expect_true(abs(mean(Q) / (p1 * p2) - 1) < 0.25, label = f)
  }
})

test_that("ind.ictest detects strong dependence", {
  set.seed(14)
  Z <- cbind(rt(150, 5), runif(150, -1, 1), rnorm(150))
  W <- cbind(Z, Z[, 1] + 0.3 * rnorm(150), rnorm(150))
  expect_lt(ind.ictest(W, 1:3, 4:5)$p.value, 1e-4)
})

test_that("independence test p-values are not rounded to zero", {
  set.seed(11)
  Z <- rmvnorm(400, rep(0, 2))
  # strong but not extreme dependence: the p-values are far below 1e-16
  # but still above the smallest representable number
  X <- cbind(Z, Z + rmvnorm(400, rep(0, 2)))
  expect_gt(ind.ctest(X, 1:2, 3:4)$p.value, 0)
  for (s in c("sign", "rank", "normal")) expect_gt(ind.ictest(X, 1:2, 3:4, scores = s)$p.value, 0)
})
