#'
#' Run spectral clustering on a (possibly directed) (possibly weighted) graph.
#'
#' It does
#' 1. do a pass-to-rank for a weighted graph (\eqn{PTR}, no-op for an unweighted graph),
#' 2. do a graph spectral embedding (\code{ASE} or \code{LSE}) with a _diagonal augmentation_,
#' 3. do a _dimension reduction_ (\eqn{ZG}) and merge left and right vectors (no-op for an undirected graph),
#' 4. cluster vertices (\eqn{GMM} or \eqn{Kmeans}).
#'
#' @param g a graph in \code{igraph} format
#' @param dmax maximum dimension for embedding
#' @param elb an index for elbow
#' @param embed either \code{ASE} or \code{LSE}, spectral embedding method
#' @param clustering either \code{GMM} or \code{Kmeans}, clustering method
#' @param use.ptr boolean to determine whether to perform pass-to-rank or not, default is \code{TRUE}
#'
#' @return \code{g} the largest connected component of the input graph
#' @return \code{mc} clustering output object
#' @return \code{Y} labels for the clustering
#' @references D.L. Sussman, M. Tang, D.E. Fishkind, and C.E. Priebe,
#' A consistent adjacency spectral embedding for stochastic blockmodel graphs,
#' Journal of the American Statistical Association, Vol. 107, No. 499, pp. 1119-1128, 2012.
#'
#' @examples
#' library(igraph)
#' data(g)
#' E(g)$weight <- runif(ecount(g), 1, 5) # add random edge weights
#' Y <- gmmase(g, dmax=20, use.ptr=TRUE, embed="ASE", clustering="Kmeans")
#'
#' @author Youngser Park <youngser@jhu.edu>
#' @export
#' @import igraph
#' @import mclust
#' @import fpc

gmmase <- function(g, dmax=20, elb=1, embed="ASE", clustering="GMM", use.ptr=TRUE)
{
#    suppressPackageStartupMessages({
#        library(igraph)
#        library(mclust)
#        library(fpc)
#    })

    cat("1. Finding an lcc...\n")
    # finding the largest connected component
    cl <- igraph::clusters(g)
    g <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
    summary(g)

    if (is.weighted(g) & use.ptr) {
        cat("2. Passing-to-rank...\n")
        g <- ptr(g)
        summary(g)
    }

    cat(paste0("3. Embedding the graph into dmax = ", dmax, "...\n"))
    if (embed=="ASE") {
        ase <- embed_adjacency_matrix(g,dmax,options=list(maxiter=10000))
    } else {
        ase <- embed_laplacian_matrix(g,dmax,options=list(maxiter=10000))
    }

    cat("4. Finding an elbow (dimension reduction)...")
    elb <- max(getElbows(ase$D)[elb],2)
    cat(", use dhat = ", elb,"\n")
    Xhat1 <- ase$X[,1:elb]
    Xhat2 <- ase$Y[,1:elb]
    Xhat <- cbind(Xhat1,Xhat2)

    cat("5. Clustering vertices...\n")
    if (clustering=="GMM") {
        mc <- Mclust(Xhat,2:9)
        plot(mc,what="BIC")
        print(summary(mc))
        Y <- mc$class
    } else {
        usepam <- ifelse(vcount(g)>2000, FALSE, TRUE)
        crit <- ifelse(vcount(g)>2000, "multiasw", "asw")
        mc <- pamk(Xhat,2:9,usepam=usepam, criterion=crit)
        plot(mc$crit, type="b")
        Y <- mc$pamobj$cluster
    }

    return(list(g=g,mc=mc,Y=Y))
}
