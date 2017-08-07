ptr <- function(g)
{
    suppressMessages(library(igraph))

#    if (class(g) == "igraph") {
#        m <- as.matrix(g[])
#    } else if (class(g) == "Matrix") {
#        m <- as.matrix(g)
#    } else m <- g

#    tmp <- m[m != 0]
#    nnz <- length(tmp)

    if (is.weighted(g)) {
        W <- E(g)$weight
    } else {
        W <- rep(1,ecount)
    }
    nnz <- length(W)
    rk <- rank(W)
    E(g)$weight <- rk * 2 / (nnz+1)
#    g <- graph.adjacency(m,weighted=TRUE)
    return(g)
}
