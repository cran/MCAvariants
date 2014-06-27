insertval2 <-
function (x, nmod) 
{
    xx <- matrix(0, length(x), nmod)
    for (i in 1:length(x)) xx[i, x[i]] <- 1
    list(xx = xx)
}
