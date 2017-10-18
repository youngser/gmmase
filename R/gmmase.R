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
#' @param Kmax maximum number of clusters
#' @param elb an index for elbow
#' @param embed either \code{ASE} or \code{LSE}, spectral embedding method
#' @param clustering either \code{GMM} or \code{Kmeans}, clustering method
#' @param use.ptr boolean to determine whether to perform pass-to-rank or not, default is \code{TRUE}
#' @param verbose boolean to determine whether to display an intermediate fitting progress status of \code{mclust} or not, default is \code{TRUE}
#' @param doplot boolean to determine whether to draw plots or not, default is \code{TRUE}
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

gmmase <- function(g, dmax=2, elb=1, lcc=TRUE, embed="ASE", clustering="GMM", Kmax=9, use.ptr=TRUE, verbose=TRUE, doplot=FALSE)
{
#    suppressPackageStartupMessages({
#        library(igraph)
#        library(mclust)
#        library(fpc)
#    })

    if (lcc) {
        cat("1. Finding an lcc...\n")
        # finding the largest connected component
        cl <- igraph::clusters(g)
        g <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
    }
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
    #    ase <- embed_laplacian_matrix(g,dmax,type="I-DAD",options=list(maxiter=1000))
    }

    cat("4. Finding an elbow (dimension reduction)...")
    elb <- max(getElbows(ase$D,plot=doplot)[elb],2)
    cat(", use dhat = ", elb,"\n")
    Xhat1 <- ase$X[,1:elb]
    if (!is.directed(g)) Xhat2 <- NULL else Xhat2 <- ase$Y[,1:elb]
    Xhat <- cbind(Xhat1,Xhat2)

    cat("5. Clustering vertices...")
    mc <- Y <- NULL
    if (clustering=="GMM") {
        if (length(Kmax)>1) {
#             for (i in 1:length(Kmax)) {
#                 mc[[i]] <- Mclust(Xhat, G=Kmax[i], verbose=verbose)
#                 Y[[i]] <- mc[[i]]$class
#                 if (doplot & i==1) {
# #                    plot(mc[[i]],what="BIC")
#                     summary(mc)
#                 }
#             }
            mc <- Mclust(Xhat, Kmax, verbose=verbose)
        } else {
            mc <- Mclust(Xhat,2:Kmax, verbose=verbose)
        }
        cat(", Khat = ", mc$G, "\n")
        if (doplot) plot(mc,what="BIC")
        print(summary(mc))
        Y <- mc$class
    } else {
        usepam <- ifelse(vcount(g)>2000, FALSE, TRUE)
        crit <- ifelse(vcount(g)>2000, "multiasw", "asw")
        if (length(Kmax)>1) {
#            for (i in 1:length(Kmax)) {
#                mc[[i]] <- pamk(Xhat, Kmax[[i]], usepam=usepam, criterion=crit)
#                Y[[i]] <- mc[[i]]$pamobj$cluster
#                if (doplot & i==1) {
##                    plot(mc[[i]]$crit, type="b")
#                    print(table(Y))
#                }
#            }
            mc <- pamk(Xhat, Kmax, usepam=usepam, criterion=crit)
        } else {
            mc <- pamk(Xhat,2:Kmax,usepam=usepam, criterion=crit)
        }
        if (doplot) plot(mc$crit, type="b")
        Y <- mc$pamobj$cluster
        cat(", Khat = ", max(Y), "\n")
        print(table(Y))
        mc$data <- Xhat
        mc$class <- Y
    }

    return(list(g=g,ase=ase,elb=elb,mc=mc,Y=Y))
}
