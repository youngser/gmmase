read.user.graph <- function(fname)
{
    dat <- as.matrix(read.table(fname, header=F, stringsAsFactors=FALSE))

    if (ncol(dat)==2) {
        cat("building a graph using an edgelist...\n")
        g <- graph.edgelist(dat, directed=TRUE)
    } else {
        cat("building a graph using an adjacency matrix...\n")
        g <- graph.adjacency(dat, directed=TRUE)
    }
    return(g)
}
