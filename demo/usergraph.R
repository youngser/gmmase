suppressPackageStartupMessages({
    library(igraph)
    library(mclust)
    library(gmmase)
})

## Please set this TRUE to use user's own graph!
user <- FALSE
if (user) {
    fname <- readline(prompt="Enter a file name (e.g., /path/edgelist.txt): ")
    g <- read_graph(fname, format="edgelist") # please read igraph manual page for details, e.g., for other types of graph format it can handle.
} else {
    data(g10000)
}
summary(g)

res <- gmmase(g)
Y <- res$class
