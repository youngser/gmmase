suppressPackageStartupMessages({
    library(igraph)
    library(mclust)
    library(gmmase)
})

print("1. Read an user graph...")

## Please uncomment these two lines and comment the line "data(g10000)") to read use graph
#fname <- readline(prompt="Enter a file name (e.g., /path/edgelist.txt): ")
#g <- read_graph(fname, format="edgelist") # please read igraph manual page for details.
data(g10000)
summary(g)

# finding the largest connected component
cl <- igraph::clusters(g)
g <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
# add edge weights
E(g)$weight <- runif(ecount(g), 1, 5)
summary(g)

print("2. Passing-to-rank...")
g <- ptr(g)
summary(g)

print("3. Embedding the graph into dmax=20...")
dmax <- 20
ase <- embed_adjacency_matrix(g,dmax,options=list(maxiter=10000))

print("4. Finding an elbow (dimension reduction)...")
elb <- max(getElbows(ase$D)[1],2)
cat("use dhat = ", elb,"\n")
Xhat1 <- ase$X[,1:elb]
Xhat2 <- ase$Y[,1:elb]
Xhat <- cbind(Xhat1,Xhat2)
dim(Xhat)

print("5. Clustering using GMM...")
mc <- Mclust(Xhat)
plot(mc,what="BIC")
summary(mc)
