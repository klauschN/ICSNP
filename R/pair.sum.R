`pair.sum` <-
function(X)
    {
     .Call(C_pairsum, num.matrix(X))
    }
