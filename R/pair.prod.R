`pair.prod` <-
function(X)
    {
     .Call(C_pairprod, num.matrix(X))
    }
