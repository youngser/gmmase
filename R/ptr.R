#'
#' Run pass-to-rank on a weighted graph.
#'
#' It extract edge weight \eqn{W} and divide \eqn{2*R / (|E|+1)} where \eqn{R} is the rank of \eqn{W} and \eqn{|E|} is the number of edges.
#'
#' @param g a graph in \code{igraph} format
#'
#' @author Youngser Park <youngser@jhu.edu>
#' @export
#'

ptr <- function(g)
{
    suppressMessages(library(igraph))

    if (is.weighted(g)) {
        W <- E(g)$weight
    } else {
        W <- rep(1,ecount(g))
    }

    E(g)$weight <- rank(W)*2 / (ecount(g)+1)
    return(g)
}
