suppressPackageStartupMessages({
    library(gmmase)
    library(igraph)
})

set.seed(123)
pm <- cbind( c(.2, .001), c(.001, .3) )
g <- sample_sbm(1000, pref.matrix=pm, block.sizes=c(300,700), directed=TRUE)
E(g)$weight <- stats::runif(ecount(g), 1, 5) # add random edge weights
summary(g)

Y <- gmmase(g, dmax=20, verbose=FALSE, doplot=TRUE)
table(Y$Y)
