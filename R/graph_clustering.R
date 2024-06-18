library("igraph")
library("Matrix")

INFTY <- 100000

#' MultipartiteClustering
#'
#' Iterative procedure to find all clusters
#' @param g The input graph
#' @return The list of clusters and their values
#' @export
MultipartiteClustering <- function(g) {

    cat("--- Multipartite graph clustering ---\n")
    global_l <- list()
    g <- GraphInit(g)
    i <- 1
    while (length(V(g)) > 0) {
        cat(paste("iteration", i, "- number of vertices:", length(V(g)), "\n"))
        g <- LinkageFunctionInit(g)
        if (AllZero(g)) {
            break
        }
        val <- MultipartiteBestCluster(g)
        l <- val[["cluster"]]
        g <- val[["graph"]]
        global_l[[i]] <- l
        g <- DeleteNodes(g, GetNodesInH(g))
        i <- i + 1
    }
    cat(paste("--- Finished:", length(global_l), "clusters found ---\n"))
    
    return(global_l)
}

#' GraphInit
#'
#' @param g The input graph
#' @return The same graph with "id" field
#' @export
GraphInit <- function(g) {

    V(g)$id <- 1:length(V(g))
#    g <- LinkageFunctionInit(g)
    return(g)
}

#' MultipartiteBestCluster
#'
#' Algorithm to find the best cluster from the current node selection
#' @param g The input graph
#' @return List containing the graph and the cluster
#' @export
MultipartiteBestCluster <- function(g) {

    Gamma <- c()
    FGamma <- -INFTY
    cat("Compute the best cluster\n")
    while (length(GetNodesInH(g)) > 0) {
        Mt <- FindArgMinLinkage(g)
        FHt <- LinkageVal(g, Mt)[1]
        all_selected <- sum(GetNodesInH(g) %in% Mt) == length(GetNodesInH(g))
        if (all_selected || AllZero(g)) {
            break
        } else {
            g <- LinkageFunctionUpdate(g = g, node.id = Mt)
            if (FHt > FGamma) {
                Gamma <- GetNodesInH(g)
                FGamma <- FHt
            }
        }
    }
    g <- SetNodesInH(g, node.id = Gamma)
    cat(paste("Solution:", length(Gamma), "nodes, value:", FGamma, "\n"))
    
    l <- list("solution" = Gamma, "value" = FGamma)
    return(list("graph" = g, "cluster" = l)) 
}

#' AllZero
#'
#' @param g The input graph
#' @return boolean value indicating whether all active nodes are 0 or not
#' @export
AllZero <- function(g) {

    zero_val <- which(LinkageVal(g, GetNodesInH(g)) == 0)
    return(length(zero_val) == length(GetNodesInH(g)))
}

#' LinkageFunctionInit
#'
#' Computes the values of the linkage function on each node of the graph
#' @param g The input graph
#' @return The same graph with "linkage.value" field
#' @export
LinkageFunctionInit <- function(g) {

    cat("Initialise the linkage function\n")
    # initialise the optimal set
    V(g)$is.inH <- rep(TRUE, length(V(g)))
    # initialise the linkage value = sum of the weights of incident edges, for every node
    V(g)$linkage.value <- sapply(V(g), function(x) sum(edge_attr(g, "weight", incident(g,x))))
    
    return(g)
}

#' LinkageFunctionUpdate
#'
#' Remove the node indexes (1-based) and update the linkage function
#' @param g The input graph
#' @param node.id node IDs to be removed from the active nodes
#' @return The same graph with "linkage.value" field updated
#' @export
LinkageFunctionUpdate <- function(g, node.id) {

    node.idx <- Id2Idx(g, node.id)
    V(g)$is.inH[node.idx] <- FALSE
    edges1 <- unlist(incident_edges(g, node.idx)) # edges incident in node.idx
    edges2 <- unlist(incident_edges(g, V(g)[V(g)$is.inH])) # edges incident in any node of H
    edges <- edges1[edges1 %in% edges2] # edges connecting node.idx and all the nodes in H
    if (length(edges) > 0) { # the graph is not complete because edges with weight = 0 are removed by LoadGraph()
        nodes <- c(ends(g, edges))
        nodes <- nodes[!(nodes %in% node.idx)] # nodes adjacent to node.idx, in the same order as the edges, with repeats
        val_edges <- edge_attr(g, "weight", edges) # weight of the edges incident to the nodes above
        idx <- unique(nodes)
        val <- sapply(idx, function(x) sum(val_edges[nodes == x])) # each entry is the sum for the nodes in H
        V(g)$linkage.value[idx] <- V(g)$linkage.value[idx] - 2*val
    }
    
    return(g)
}

#' FindArgMinLinkage
#'
#' Find the node with the minimum linkage value
#' @param g The input graph
#' @return ID of the node with min linkage value
#' @export
FindArgMinLinkage <- function(g) {

    val <- min(V(g)$linkage.value[V(g)$is.inH])
    idx <- as_ids(V(g)[V(g)$linkage.value <= val & V(g)$is.inH]) 
    return(Idx2Id(g, idx))
}

#' LinkageVal
#'
#' @param g The input graph
#' @param node.id the ID of the node(s)
#' @return linkage value for node.id
#' @export
LinkageVal <- function(g, node.id) {

    node.idx <- Id2Idx(g, node.id)
    v <- c(node.idx)
    return(V(g)$linkage.value[node.idx])
}

#' DeleteNodes
#'
#' @param g The input graph
#' @param node.id the ID of the node(s)
#' @return graph with node.id deleted
#' @export
DeleteNodes <- function(g, node.id) {

    node.idx <- Id2Idx(g, node.id)
    g <- delete_vertices(g, node.idx)
    return(g)
}

#' GetNodesInH
#'
#' @param g The input graph
#' @return The sequence of active nodes
#' @export
GetNodesInH <- function(g) {

    idx <- as_ids(V(g)[V(g)$is.inH])
    return(Idx2Id(g, idx))
}

#' GetNodesInH
#'
#' @param g The input graph
#' @param node.id The sequence of nodes
#' @return the same graph with node.id nodes active
#' @export
SetNodesInH <- function(g, node.id) {

    idx <- Id2Idx(g, node.id)
    V(g)$is.inH <- as_ids(V(g)) %in% idx
    return(g)
}

#' Id2Idx
#'
#' @param g The input graph
#' @param node.id The sequence of node IDs
#' @return node.idx corresponding to node.id
#' @export
Id2Idx <- function(g, node.id) {

    if (sum(V(g)$id %in% node.id) < length(node.id)) {
        stop(paste("Node id", node.id, "not found in this graph\n"))
    }
    return(as_ids(V(g)[match(node.id,V(g)$id)]))
}

#' IdxId
#'
#' @param g The input graph
#' @param node.idx The sequence of node indexes (1-based)
#' @return node.id corresponding to node.idx
#' @export
Idx2Id <- function(g, node.idx) {

    if (length(node.idx) == 0) {
        return(NULL)
    }
    if (max(node.idx) > length(V(g))) {
        stop(paste("Node idx", node.idx, "exceeds graph size\n"))
    }
    return(V(g)$id[node.idx])
}


