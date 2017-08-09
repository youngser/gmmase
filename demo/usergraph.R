suppressPackageStartupMessages({
    library(igraph)
    library(fpc)
    library(gmmase)
})

## Please set this TRUE to use user's own graph!
user <- FALSE
if (user) {
#    fname <- readline(prompt="Enter a file name (e.g., /path/edgelist.txt): ")
    fname <- "g-10000.el"
    g <- read_graph(fname, format="edgelist")
} else {
    data(g10000)
}

E(g)$weight <- runif(ecount(g), 1, 5) # add random edge weights
summary(g)

# GMM for a large graph takes long time (e.g., it takes >5 minutes for a graph with 10K vertices), so "kmeans" is recommended.
system.time(Y <- gmmase(g, embed = "LSE", clustering="kmeans"))
table(Y)
