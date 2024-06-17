
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <file_list> <out>\n")
    q()
}

FILE_LIST <- args[1]
OUT_PREFIX <- args[2]

library("multiGraphClust")

file.list <- as.list(read.table(FILE_LIST)[,1])
val <- ParseMultipartiteGraphList(file.list)
M <- val[["matrix"]]
g <- LoadGraph(M)

block_id <- val[["block_id"]]
out <- paste0(OUT_PREFIX, "_block_id.txt")
write(block_id, file = out)

# cluster using the linkage function (Vashist et al, 2007)
# here, the distance between node groups (p() function) is assumed to be constant = 1
global_l <- MultipartiteClustering(g)

# save to file
out <- paste0(OUT_PREFIX, "_out.txt")
PrintMultipartiteClustering(cluster.list = global_l, file.name = out)

q()

