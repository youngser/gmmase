gmmase <- function(g, dmax=20, embed="ASE", clustering="GMM", add.weight=FALSE, use.ptr=TRUE)
{
    suppressPackageStartupMessages({
        library(igraph)
        library(mclust)
        library(fpc)
    })

    cat("1. Finding an lcc...\n")
    # finding the largest connected component
    cl <- igraph::clusters(g)
    g <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))

    if (add.weight) {
        # add edge weights
        E(g)$weight <- runif(ecount(g), 1, 5)
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
    }

    cat("4. Finding an elbow (dimension reduction)...")
    elb <- max(getElbows(ase$D)[1],2)
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

    return(Y)
}
