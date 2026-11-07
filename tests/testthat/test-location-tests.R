test_that("HotellingsT2 one-sample test agrees with the textbook formula", {
  data("pulmonary", package = "ICSNP", envir = environment())
  X <- as.matrix(pulmonary)
  n <- nrow(X); p <- ncol(X)
  xbar <- colMeans(X)
  T2 <- n * drop(t(xbar) %*% solve(cov(X)) %*% xbar)
  Fstat <- T2 * (n - p) / (p * (n - 1))

  res.f <- HotellingsT2(X)
  expect_s3_class(res.f, "htest")
  expect_equal(unname(res.f$statistic), Fstat, tolerance = 1e-10)
  expect_equal(unname(res.f$parameter), c(p, n - p))
  expect_equal(res.f$p.value, pf(Fstat, p, n - p, lower.tail = FALSE), tolerance = 1e-10)

  res.chi <- HotellingsT2(X, test = "chi")
  expect_equal(unname(res.chi$statistic), T2, tolerance = 1e-10)
  expect_equal(res.chi$p.value, pchisq(T2, p, lower.tail = FALSE), tolerance = 1e-10)
})

test_that("HotellingsT2 two-sample test agrees with the textbook formula and the formula method", {
  set.seed(1)
  X <- rmvnorm(20, c(0, 0, 0))
  Y <- rmvnorm(25, c(0.5, 0, 0))
  n1 <- nrow(X); n2 <- nrow(Y); p <- ncol(X)
  S <- ((n1 - 1) * cov(X) + (n2 - 1) * cov(Y)) / (n1 + n2 - 2)
  d <- colMeans(X) - colMeans(Y)
  T2 <- n1 * n2 / (n1 + n2) * drop(t(d) %*% solve(S) %*% d)
  Fstat <- T2 * (n1 + n2 - p - 1) / (p * (n1 + n2 - 2))

  res <- HotellingsT2(X, Y)
  expect_equal(unname(res$statistic), Fstat, tolerance = 1e-10)
  expect_null(dim(res$statistic))
  expect_named(res$statistic, "T.2")
  expect_null(dim(HotellingsT2(X, Y, test = "chi")$statistic))
  expect_equal(unname(res$parameter), c(p, n1 + n2 - p - 1))

  Z <- rbind(X, Y)
  g <- factor(rep(1:2, c(n1, n2)))
  expect_equal(HotellingsT2(Z ~ g)$statistic, res$statistic)
})

test_that("HotellingsT2 gives small p-values without rounding them to zero", {
  set.seed(9)
  X <- rmvnorm(50, c(3, 3))
  Y <- rmvnorm(50, c(0, 0))
  expect_gt(HotellingsT2(X)$p.value, 0)
  expect_gt(HotellingsT2(X, Y)$p.value, 0)
  expect_gt(HotellingsT2(X, Y, test = "chi")$p.value, 0)
})

test_that("formula methods use the data, subset and na.action arguments", {
  ir <- iris[iris$Species != "virginica", ]
  X <- as.matrix(ir[ir$Species == "setosa", 1:2])
  Y <- as.matrix(ir[ir$Species == "versicolor", 1:2])
  fo <- cbind(Sepal.Length, Sepal.Width) ~ Species

  expect_equal(HotellingsT2(fo, data = ir)$statistic, HotellingsT2(X, Y)$statistic)
  expect_equal(HotellingsT2(fo, data = iris, subset = Species != "virginica")$statistic,
               HotellingsT2(X, Y)$statistic)
  expect_equal(rank.ctest(fo, data = ir)$statistic, rank.ctest(X, Y)$statistic)
  expect_equal(rank.ctest(fo, data = ir, scores = "sign")$statistic,
               rank.ctest(X, Y, scores = "sign")$statistic)
  # c samples
  Z3 <- as.matrix(iris[, 1:2])
  expect_equal(rank.ctest(fo, data = iris)$statistic, rank.ctest(Z3 ~ iris$Species)$statistic)
  # a numeric matrix as data
  M <- cbind(a = c(X[, 1], Y[, 1]), b = c(X[, 2], Y[, 2]), g = rep(1:2, c(50, 50)))
  expect_equal(HotellingsT2(cbind(a, b) ~ g, data = M)$statistic, HotellingsT2(X, Y)$statistic)

  # missing values: the default is to fail, as for the default methods
  irNA <- ir
  irNA$Sepal.Length[3] <- NA
  expect_error(HotellingsT2(fo, data = irNA), "missing values")
  expect_error(rank.ctest(fo, data = irNA), "missing values")
  expect_equal(HotellingsT2(fo, data = irNA, na.action = na.omit)$statistic,
               HotellingsT2(X[-3, ], Y)$statistic)
  expect_equal(rank.ctest(fo, data = irNA, na.action = na.omit)$statistic,
               rank.ctest(X[-3, ], Y)$statistic)
})

test_that("HotellingsT2 and HP.loc.test are affine invariant", {
  set.seed(2)
  X <- rmvnorm(40, c(0.3, 0, 0), diag(c(1, 2, 3)))
  XA <- X %*% t(A3)
  expect_equal(HotellingsT2(XA)$statistic, HotellingsT2(X)$statistic, tolerance = 1e-8)
  for (s in c("sign", "rank", "normal")) {
    expect_equal(HP.loc.test(XA, score = s)$statistic,
                 HP.loc.test(X, score = s)$statistic, tolerance = 1e-4)
  }
})

test_that("HP.loc.test returns chi-square p-values and valid permutation p-values", {
  set.seed(3)
  X <- rmvnorm(30, c(0, 0, 0))
  res <- HP.loc.test(X)
  expect_s3_class(res, "htest")
  expect_equal(res$p.value, pchisq(unname(res$statistic), 3, lower.tail = FALSE))
  res.perm <- HP.loc.test(X, method = "permutation", n.perm = 50)
  expect_true(res.perm$p.value >= 0 && res.perm$p.value <= 1)
  expect_equal(unname(res.perm$parameter), 50)
})

test_that("univariate marginal signed-rank test agrees with the textbook formula", {
  set.seed(4)
  x <- rnorm(25, 0.3)
  r <- rank(abs(x))
  Tstat <- sum(sign(x) * r)^2 / sum(r^2)
  res <- rank.ctest(matrix(x), scores = "rank")
  expect_equal(unname(res$statistic), Tstat, tolerance = 1e-10)
  expect_equal(res$p.value, pchisq(Tstat, 1, lower.tail = FALSE), tolerance = 1e-10)
})

test_that("marginal tests are invariant under coordinatewise rescaling, sign changes and permutations", {
  set.seed(5)
  X <- cbind(rt(30, 3), rnorm(30, 0.4), runif(30, -1, 2))
  D <- diag(c(2, -0.5, 10))
  for (s in c("sign", "rank", "normal")) {
    ref <- rank.ctest(X, scores = s)$statistic
    expect_equal(rank.ctest(X %*% D, scores = s)$statistic, ref, tolerance = 1e-8)
    expect_equal(rank.ctest(X[, 3:1], scores = s)$statistic, ref, tolerance = 1e-8)
  }
})

test_that("rank.ctest handles one, two and c samples consistently", {
  set.seed(6)
  X <- rmvnorm(20, c(0, 0))
  Y <- rmvnorm(25, c(1, 0))
  W <- rmvnorm(15, c(0, 1))
  Z <- rbind(X, Y)
  g <- factor(rep(1:2, c(20, 25)))
  for (s in c("sign", "rank", "normal")) {
    res2 <- rank.ctest(X, Y, scores = s)
    expect_equal(unname(res2$parameter), 2)
    expect_equal(rank.ctest(Z ~ g, scores = s)$statistic, res2$statistic)
    # the c-sample code with two groups gives the two-sample statistic
    expect_equal(ICSNP:::MaRaTe.internal.csample(Z, g, scores = s)$test.statistic,
                 unname(res2$statistic), tolerance = 1e-8)
  }
  Z3 <- rbind(X, Y, W)
  g3 <- factor(rep(1:3, c(20, 25, 15)))
  res3 <- rank.ctest(Z3 ~ g3)
  expect_s3_class(res3, "htest")
  expect_equal(unname(res3$parameter), 2 * (3 - 1))
  expect_equal(res3$p.value, pchisq(unname(res3$statistic), 4, lower.tail = FALSE))
})

test_that("rank.ictest gives the same result for an ics object and its components", {
  set.seed(7)
  Z <- cbind(rt(60, 5), runif(60, -1, 1), rnorm(60))
  X <- Z %*% t(A3)
  ics.X <- ics(X, S1 = cov, S2 = cov4)
  comps <- as.matrix(ics.components(ics.X))
  for (s in c("sign", "rank", "normal")) {
    expect_equal(rank.ictest(ics.X, scores = s)$statistic,
                 rank.ictest(comps, scores = s)$statistic)
  }
  set.seed(8); p1 <- rank.ictest(comps, method = "permutation", n.simu = 100)$p.value
  set.seed(8); p2 <- rank.ictest(comps, method = "permutation", n.simu = 100)$p.value
  expect_identical(p1, p2)
  expect_true(p1 >= 0 && p1 <= 1)
})

test_that("chi-square based p-values are not rounded to zero", {
  set.seed(10)
  X <- rmvnorm(300, c(3, 3, 3))
  expect_gt(rank.ctest(X)$p.value, 0)
  expect_gt(HP.loc.test(X)$p.value, 0)
  expect_gt(rank.ictest(X)$p.value, 0)
  Y <- rmvnorm(300, c(0, 0, 0))
  expect_gt(rank.ctest(X, Y)$p.value, 0)
  Z <- rbind(X, Y, rmvnorm(300, c(-3, -3, -3)))
  g <- factor(rep(1:3, each = 300))
  expect_gt(rank.ctest(Z ~ g)$p.value, 0)
})
