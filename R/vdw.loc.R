# signed-rank statistic with van der Waerden scores,
# T(mu) = sum sign(y_i - mu) qnorm(1/2 + R_i^+(mu) / (2 (n + 1))),
# where R_i^+(mu) is the rank of |y_i - mu|. T is nonincreasing in mu and
# constant between neighbouring Walsh averages of y.
T.vdw.loc <- function(mu, y)
    {
    n <- length(y)
    z <- y - mu
    sum(sign(z) * qnorm(rank(abs(z)) / (2 * (n + 1)) + 0.5))
    }

# The estimate is the Hodges-Lehmann estimate based on T,
# (sup{mu : T(mu) > 0} + inf{mu : T(mu) < 0}) / 2.
# The arguments int.diff and maxiter are no longer used.
vdw.loc <- function(x, int.diff = 10, maxiter = 1000, na.action = na.fail)
    {
    if (!is.vector(x)) stop("'x' must be a numeric vector")
    x <- na.action(x)
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    x <- as.vector(x)
    n <- length(x)
    if (n == 0) stop("'x' must contain at least one observation")

    # the distinct Walsh averages w[1] < ... < w[M]
    w <- sort(c(x, pair.sum(as.matrix(x)) / 2))
    w <- w[c(TRUE, diff(w) > 0)]
    M <- length(w)
    if (M == 1) return(w)

    # T is constant on (w[k], w[k+1]), so it is evaluated at the midpoints
    mids <- (w[-1] + w[-M]) / 2
    # differences in T below tol are rounding errors
    tol <- 1000 * .Machine$double.eps * sum(qnorm((1:n) / (2 * (n + 1)) + 0.5))

    # bisection for the last midpoint index k in 0, ..., M - 1 with pred TRUE,
    # where pred is TRUE up to some k and FALSE after it (index 0 stands for
    # -Inf, where T is positive)
    last.true <- function(pred)
        {
        lo <- 0
        hi <- M
        while (hi - lo > 1)
            {
            mid <- (lo + hi) %/% 2
            if (pred(T.vdw.loc(mids[mid], x))) lo <- mid else hi <- mid
            }
        lo
        }
    # T > 0 on the intervals up to k.pos, so sup{T > 0} = w[k.pos + 1]
    k.pos <- last.true(function(t) t > tol)
    # T < 0 on the intervals after k.nonneg, so inf{T < 0} = w[k.nonneg + 1]
    k.nonneg <- last.true(function(t) t >= -tol)
    (w[k.pos + 1] + w[k.nonneg + 1]) / 2
    }
