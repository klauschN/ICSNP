symm.huber <- function(X, qg=0.9, init=NULL, eps = 1e-06, maxiter = 100, na.action=na.fail)
    {
    X <- na.action(X)
    X <- as.matrix(X)
    CNAMES <- colnames(X)
    p <- ncol(X)
    
    c.square <- 2 * qchisq(qg, p)
    sigma.square <- 2 * pchisq(c.square / 2, p + 2) + (c.square / p) * (1 - qg)
    
    if (p<2) stop("'X' must be at least bivariate")  

    n<-dim(X)[1]

    if (is.numeric(init)) V.0 <- init else V.0<- cov(X)

    # the fixed-point iteration is done in C (see src/icsnp.c)
    V.new <- .Call(C_symm_huber, num.matrix(X), num.matrix(V.0), as.double(c.square),
                   as.double(sigma.square), as.double(eps), as.double(maxiter))
    attr(V.new, "iter") <- NULL
    colnames(V.new) <- CNAMES
    rownames(V.new) <- CNAMES
    V.new
    }

# Internal function (used in the tests): one step of the iteration, the
# weighted average over all pairs of the outer products of the pairwise
# differences, where V is the current scatter estimate
SSCov.hub<-function(X,V,c.s,sig.s)
{
.Call(C_symm_huber_step, num.matrix(X), num.matrix(V), as.double(c.s), as.double(sig.s))
}
