################################################################################
## Code reproducing the illustrative examples of
##   Nordhausen, K., Oja, H. and Tyler, D. E.
##   "ICSNP: An R package for multivariate nonparametric location tests and
##    robust estimates of location and shape" (submitted).
##
## The script reproduces all results of Section 3 (Illustrative examples) in
## the order in which they appear in the paper, and it creates Figure 1 as the
## file fig_shapes.pdf in the working directory. The results in the paper were
## obtained with R 4.7.0 (development version) and ICSNP 1.2-0 (see the output
## of sessionInfo() at the end). Only the computation times in Example 5 are
## expected to differ between computers.
################################################################################

library("ICSNP")     # also loads the packages mvtnorm and ICS
options(digits = 4)  # the printed test output in the paper uses 4 digits
cat("ICSNP version:", as.character(packageVersion("ICSNP")), "\n")


## =============================================================================
## Example 1: one-sample location tests for the pulmonary data
## =============================================================================

## Changes in three measurements of the pulmonary function of 12 workers after
## six hours of exposure to cotton dust. The null hypothesis is that the
## location of the changes is the origin.
data("pulmonary")

## marginal signed-rank test (the output shown in the paper)
rank.ctest(pulmonary, scores = "rank")

## Hotelling's T^2 test, the marginal sign and normal score tests, and the
## affine invariant sign, signed-rank and van der Waerden tests of Hallin and
## Paindaveine
tests.pulm <- list(
  "Hotelling's T2"        = HotellingsT2(pulmonary),
  "marginal sign"         = rank.ctest(pulmonary, scores = "sign"),
  "marginal signed-rank"  = rank.ctest(pulmonary, scores = "rank"),
  "marginal normal score" = rank.ctest(pulmonary, scores = "normal"),
  "HP sign"               = HP.loc.test(pulmonary, score = "sign"),
  "HP signed-rank"        = HP.loc.test(pulmonary, score = "rank"),
  "HP van der Waerden"    = HP.loc.test(pulmonary, score = "normal"))

## the p-values reported in the paper
cat("\np-values of the one-sample tests:\n")
print(round(sapply(tests.pulm, function(x) x$p.value), 3))

## Affine invariance: transform the observations as x_i -> A x_i with a
## nonsingular matrix A. The statistics of the invariant tests do not change,
## the marginal signed-rank statistic does.
A <- matrix(c(1, 2, 0, 0.5, 1, 3, 1, -1, 2), ncol = 3)
pulm.A <- as.matrix(pulmonary) %*% t(A)
cat("\nTest statistics for the original and the transformed data:\n")
stat.A <- rbind(
  "Hotelling's T2"       = c(HotellingsT2(pulmonary)$statistic,
                             HotellingsT2(pulm.A)$statistic),
  "HP signed-rank"       = c(HP.loc.test(pulmonary, score = "rank")$statistic,
                             HP.loc.test(pulm.A, score = "rank")$statistic),
  "marginal signed-rank" = c(rank.ctest(pulmonary, scores = "rank")$statistic,
                             rank.ctest(pulm.A, scores = "rank")$statistic))
colnames(stat.A) <- c("original", "transformed")
print(round(stat.A, 3))


## =============================================================================
## Example 2: c-sample location test for the iris data
## =============================================================================

## marginal rank test for the equality of the locations of the three species,
## using the formula interface
rank.ctest(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species,
           data = iris, scores = "rank")


## =============================================================================
## Example 3: robust estimates of location and shape (Figure 1)
## =============================================================================

## 95 observations from a bivariate normal distribution with correlation 0.7
## and five outliers around (6, -6)
set.seed(123)
Sigma <- matrix(c(1, 0.7, 0.7, 1), 2)
X <- rbind(rmvnorm(95, c(0, 0), Sigma),
           rmvnorm(5, c(6, -6), diag(0.1, 2)))
colnames(X) <- c("X1", "X2")

HR <- HR.Mest(X)            # Hettmansperger-Randles location and shape
D  <- duembgen.shape(X)     # Duembgen's shape matrix
SH <- symm.huber(X)         # symmetrized Huber scatter matrix
m  <- spatial.median(X)     # spatial median

## To compare only the shapes, all scatter matrices are standardized to have
## determinant one (duembgen.shape is already standardized in this way).
shp <- function(S) S / det(S)^(1 / ncol(S))
est <- list(
  "Sample covariance"      = list(center = colMeans(X), shape = shp(cov(X))),
  "Hettmansperger-Randles" = list(center = HR$center,   shape = shp(HR$scatter)),
  "Dümbgen"           = list(center = m,           shape = D),
  "Symmetrized Huber"      = list(center = m,           shape = shp(SH)))

for (nm in names(est)) {
  cat("\n", nm, ": center =", round(est[[nm]]$center, 3), "\n")
  print(round(est[[nm]]$shape, 3))
}
cat("\nTrue shape:\n"); print(round(shp(Sigma), 3))

## the correlation implied by the sample covariance matrix (negative because
## of the outliers) and by the robust estimates
cat("\nImplied correlations:\n")
print(round(sapply(est, function(e) cov2cor(e$shape)[1, 2]), 2))

## Figure 1: the data and the contour ellipses
## {x : (x - mu)' V^{-1} (x - mu) = qchisq(0.95, 2)} for the estimates mu and V.
## Each shape matrix V is rescaled to the trace of the true covariance matrix
## (trace 2), so that only the location and the shape are compared.
ellipse <- function(center, S, n.points = 200) {
  t <- seq(0, 2 * pi, length.out = n.points)
  S <- S * qchisq(0.95, 2) / mean(diag(S)) * mean(diag(Sigma))
  t(center + t(cbind(cos(t), sin(t)) %*% chol(S)))
}
## Okabe-Ito colors, which are distinguishable for colorblind readers; the line
## types keep the ellipses distinguishable also in grayscale
cols <- c("#D55E00", "#0072B2", "#009E73", "#CC79A7")
ltys <- c(2, 1, 3, 4)

pdf("fig_shapes.pdf", width = 5.5, height = 5.5)
par(mar = c(4, 4, 1, 1))
plot(X, pch = 19, col = rep(c("grey65", "black"), c(95, 5)), cex = 0.7,
     xlim = c(-4, 7), ylim = c(-7, 4), las = 1)
for (i in seq_along(est))
  lines(ellipse(est[[i]]$center, est[[i]]$shape), lty = ltys[i], col = cols[i], lwd = 3)
legend("bottomleft", names(est), lty = ltys, col = cols, lwd = 3, seg.len = 3,
       bty = "n", cex = 0.85)
dev.off()


## =============================================================================
## Example 4: location test in the independent component model
## =============================================================================

## Three independent components with t_5, uniform and symmetrized exponential
## distributions (the latter two standardized to unit variance), a location
## shift of 0.2 in the first component, mixed with a nonsingular matrix Omega
set.seed(2026)
n <- 200
Z <- cbind(rt(n, 5),
           runif(n, -sqrt(3), sqrt(3)),
           rexp(n) * sample(c(-1, 1), n, replace = TRUE))
Omega <- matrix(c(3, 2, 1, 2, 4, -0.5, 1, -0.5, 2), ncol = 3)
Xic <- (Z + matrix(c(0.2, 0, 0), n, 3, byrow = TRUE)) %*% t(Omega)

## invariant coordinates based on the covariance matrix and the matrix of
## fourth moments, followed by the signed-rank test of Nordhausen, Oja and
## Paindaveine (2009) and the normal score test with permutation p-value
ics.X <- ics(Xic, S1 = cov, S2 = cov4)
rank.ictest(ics.X, scores = "rank")
rank.ictest(ics.X, scores = "normal", method = "permutation", n.simu = 2000)


## =============================================================================
## Example 5: computation times of the C and R implementations of Duembgen's
## shape matrix
## =============================================================================

## n = 1000 observations from a five-variate normal distribution, that is,
## choose(1000, 2) = 499500 pairwise differences. The times depend on the
## computer used.
set.seed(1)
Xbig <- rmvnorm(1000, rep(0, 5))
print(system.time(D.C <- duembgen.shape(Xbig)))               # C implementation
print(system.time(D.R <- duembgen.shape(Xbig, in.R = TRUE)))  # pure R version
cat("largest absolute difference between the two estimates:",
    format(max(abs(D.C - D.R)), digits = 2), "\n")


## =============================================================================
## Computing environment
## =============================================================================
sessionInfo()
