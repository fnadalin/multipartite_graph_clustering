
library("multiGraphClust")
setwd("example")

file.list <- as.list(read.table("matrix_list.txt")[,1])
val <- ParseMultipartiteGraphList(file.list)
M <- val[["matrix"]]
g <- LoadGraph(M)
block_id <- val[["block_id"]]
write(block_id, file = "block_id.txt")

# cluster using the linkage function (Vashist et al, 2007)
# here, the distance between node groups (p() function) is assumed to be constant = 1

global_l <- MultipartiteClustering(g)
PrintMultipartiteClustering(cluster.list = global_l, file.name = "out.txt")

# compare to general clustering algorithms

c <- cluster_optimal(g)
m <- modularity(g, membership = membership(c))

c <- cluster_louvain(g, resolution = 1.5)
m <- modularity(g, membership = membership(c))


