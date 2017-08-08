suppressPackageStartupMessages({
    library(igraph)
    library(mclust)
    library(gmmase)
})

## Please set this TRUE to use user's own graph!
user <- TRUE
if (user) {
#    fname <- readline(prompt="Enter a file name (e.g., /path/edgelist.txt): ")
    fname <- "~/Dropbox/gmmase/data/g-10000.el"
    g <- read_graph(fname, format="edgelist")
} else {
    data(g10000)
}
summary(g)

system.time(res <- gmmase(g, add.weight=TRUE))
Y <- res$class
