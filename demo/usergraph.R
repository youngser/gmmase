suppressPackageStartupMessages({
    library(igraph)
    library(mclust)
    library(gmmase)
})

## Please uncomment these two lines and comment the line "data(g10000)") to read use graph
#fname <- readline(prompt="Enter a file name (e.g., /path/edgelist.txt): ")
#g <- read_graph(fname, format="edgelist") # please read igraph manual page for details, e.g., for other types of graph format it can handle.
data(g10000)
summary(g)

res <- gmmase(g, dmax=2)
Y <- res$class
