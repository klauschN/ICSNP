`duembgen.shape` <-
function(X,init=NULL,steps=Inf,eps=1e-6,maxiter=100,in.R=FALSE,na.action=na.fail,...)
    {
    X<-na.action(X)
    X<-as.matrix(X)
    
    n <- dim(X)[1] 
    p <- dim(X)[2]
    if (p<2) stop("'X' must be at least bivariate")  
    
    if (in.R==FALSE){
            # zero pairwise differences (tied observations) are skipped in the C code
            if (anyDuplicated(X) > 0)
                {
                id <- match(data.frame(t(X)), data.frame(t(unique(X))))
                n.del <- sum(choose(tabulate(id), 2))
                warning(paste(n.del, "pairwise difference(s) equal to zero were removed"))
                }
            if (is.null(init)) init <- cov(X)
            init <- shape.det(init)
            if(is.finite(steps)) maxiter <- Inf
            # the fixed-point iteration is done in C (see src/icsnp.c); with
            # steps = 0 the (standardized) initial value is returned
            duembgen <- .Call(C_duembgen_shape, num.matrix(X), num.matrix(init),
                              as.double(eps), as.double(maxiter), as.double(steps))
            attr(duembgen, "iter") <- NULL
            colnames(duembgen) <- colnames(X)
            rownames(duembgen) <- colnames(X)
            }
            
    
    if (in.R==TRUE){
        data2<-pair.diff(X)
        colnames(data2) <- colnames(X)
        duembgen<-tyler.shape(data2,0,init,steps,eps,maxiter,in.R,...)}
        
    return(duembgen)
    }


# Function to standardize a scatter matrix to have det 1, internal for duembgen.shape
shape.det <- function(M) M/det(M)^(1/dim(M)[2])

# Internal function (used in the tests): the average over all pairs of the
# outer products of the spatial signs of the pairwise differences. Zero
# differences are skipped, and their number is in the attribute "n.zero".
SSCov <- function(X) .Call(C_sum_diff_sign_outers, num.matrix(X))
