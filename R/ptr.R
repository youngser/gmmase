ptr <- function(g)
{
    suppressMessages(library(igraph))

    if (class(g) == "igraph") {
        m <- as.matrix(g[])
    } else if (class(g) == "Matrix") {
        m <- as.matrix(m)
    } else m <- g

    tmp <- m[m != 0]
    nnz <- length(tmp)
    rk <- rank(tmp)
    m[m != 0] <- rk * 2 / (nnz+1)
    g <- graph.adjacency(m,weighted=TRUE)
    return(g)
}
