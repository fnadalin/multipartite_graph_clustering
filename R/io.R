library("tools")
library("igraph")
library("Matrix")

EPSILON <- 0.0000001

#' ParseMultipartiteGraphList
#'
#' Parse a named, ordered list of bipartite graph file names and store it into a 0-diagonal, symmetric block matrix
#' The matrices must be parsed in row order and then column order
#' N.B. 1) the # of columns of vertex set i must equate the # of rows of vertex set i+1
#'      2) rownames, colnames: auto-detected, but not used
#' @param file.list txt file containing the matrices
#' @return the upper-triangular block matrix and the block ID for each row/column
#' @export
ParseMultipartiteGraphList <- function(file.list) {

    # get the nunber of blocks per row
    tot.blocks <- length(file.list)
    row.blocks <- (1 + sqrt(1 + 8*tot.blocks)) / 2
    if (row.blocks > floor(row.blocks) + EPSILON) {
        stop("The number of blocks is not compatible with a triangular block matrix\n")
    }

    # get the matrix dimensions
    block_id <- NULL
    i <- 1
    while (i < row.blocks) {
        f <- file.list[[i]]
        is.sparse <- file_ext(f) == "mtx"
        if (is.sparse) {
            data <- readMM(f)
        } else {
            data <- read.table(file = f, sep = "\t")
        }
        i <- i + 1
        block_id <- c(block_id, rep(i,ncol(data)))
    }
    block_id <- c(rep(1,nrow(data)), block_id)
    tot <- length(block_id)
    
    # fill the upper triangular block matrix
    M <- matrix(0, nrow = tot, ncol = tot)
    start.row <- start.col <- start <- 0
    i <- 1
    for (n in (row.blocks-1):1) {
        prev.col <- prev.col.tmp <- NA
        for (b in 1:n) {
            f <- file.list[[i]]
            is.sparse <- file_ext(f) == "mtx"
            if (is.sparse) {
                data <- readMM(f)
            } else {
                data <- as.matrix(read.table(file = f, sep = "\t"))
            }
            if (!is.na(prev.col) & prev.col != nrow(data)) { # check that the dimensions of the blocks are consistent
                stop("Column number of previous matrix and row number of current matrix differ\n")
            }
            if (b == 1) {
                start.col <- start + nrow(data)
                prev.col.tmp <- ncol(data)
            } 
            M[start.row + 1:nrow(data), start.col + 1:ncol(data)] <- data
            start.col <- start.col + ncol(data)
            i <- i + 1
        }
        prev.col <- prev.col.tmp
        start.row <- start <- start + nrow(data)
    }
    
    # symmetrise
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
    
    return(list("matrix" = M, "block_id" = block_id))
}

#' LoadGraph
#'
#' Load a matrix or a sparse matrix (Matrix object) and store it into an igraph object
#' Load also the metadata associated to the nodes (i.e. the vertex sets of the multi-partite graph)
#' @param M the upper-triangular block matrix created with ParseMultipartiteGraphList
#' @return a graph structure
#' @export
LoadGraph <- function(M) {

    # also checks that the matrix is symmetric
    g <- graph_from_adjacency_matrix(adjmatrix = M, mode = "undirected", weighted = TRUE, add.colnames = NA, add.rownames = NA)
    return(g)
}

# print the solution of the algorithm to file
# the solution can be an igraph object itself (induced by H, the output subset of nodes)
WriteGraph <- function(G, filename) {

    M <- as.matrix(g, "adjacency") # FIXME: this removes the weights!!!
    is.sparse <- file_ext(filename) == "mtx"
    if (!is.sparse) {
        write.table(x = M, file = filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
        writeMM(obj = M, file = filename)
    }
}

#' PrintMultipartiteClustering
#'
#' Print the clusters, one per line
#' @param cluster.list The output of MultipartiteClustering
#' @file.name output txt file
#' @return nothing
#' @export
PrintMultipartiteClustering <- function(cluster.list, file.name) {

    v <- sapply(cluster.list, function(x) paste(unlist(x[["solution"]]), collapse = ","))
    write(v, file = file.name)
}


