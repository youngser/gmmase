#'
#' Run pass-to-rank on a weighted graph.
#'
#' It extracts (non-zero) edge weight vector \eqn{W} from a graph and replaces it with \eqn{2*R / (|E|+1)} where \eqn{R} is the rank of \eqn{W} and \eqn{|E|} is the number of edges. This does 'no-op' for an unweighted graph.
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
