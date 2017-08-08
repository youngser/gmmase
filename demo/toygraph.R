suppressPackageStartupMessages({
    library(igraph)
    library(mclust)
    library(gmmase)
})

set.seed(123)
pm <- cbind( c(.2, .001), c(.001, .3) )
g <- sample_sbm(100, pref.matrix=pm, block.sizes=c(30,70), directed=TRUE)
summary(g)

res <- gmmase(g, add.weight=TRUE)
Y <- res$class # cluster labels
