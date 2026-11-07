`pair.diff` <-
function(X)
    {
     .Call(C_pairdiff, num.matrix(X))
    }
