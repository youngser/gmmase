suppressPackageStartupMessages({
    library(gmmase)
    library(igraph)
})

## Please set this TRUE to use user's own graph!
user <- FALSE
if (user) {
#    fname <- readline(prompt="Enter a file name (e.g., /path/edgelist.txt): ")
    fname <- "g.el"
    g <- read_graph(fname, format="edgelist")
} else {
    data(g)
}

E(g)$weight <- stats::runif(ecount(g), 1, 5) # add random edge weights
summary(g)

# GMM for a large graph takes long time (e.g., it takes >5 minutes for a graph with 10K vertices), so "kmeans" is recommended.
system.time(Y <- gmmase(g, dmax=20, embed = "LSE", clustering="kmeans", verbose=FALSE, doplot=TRUE))
table(Y$Y)
