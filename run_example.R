
setwd("example")

file.list <- as.list(read.table("matrix_list.txt")[,1])
val <- ParseMultipartiteGraphList(file.list)
M <- val[["matrix"]]
block_id <- val[["block_id"]]
g <- LoadGraph(M)

# cluster using the linkage function (Vashist et al, 2007)
# here, the distance between node groups (p() function) is assumed to be constant = 1

global_l <- MultipartiteClustering(g)

# compare to general clustering algorithms

c <- cluster_optimal(g)
m <- modularity(g, membership = membership(c))

c <- cluster_louvain(g, resolution = 1.5)
m <- modularity(g, membership = membership(c))


