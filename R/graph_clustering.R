library("igraph")
library("Matrix")

INFTY <- 100000

# iterative procedure to find all clusters
MultipartiteClustering <- function(g) {

    all_zero <- FALSE
    global_l <- list()
    g <- GraphInit(g)
    i <- 1
    while (length(V(g)) > 0 || !all_zero) {
        cat(paste("length(V(g)):", length(V(g)), "\n"))
        val <- MultipartiteBestCluster(g)
        l <- val[["cluster"]]
        g <- val[["graph"]]
        global_l[[i]] <- l
        g <- DeleteNodes(g, GetNodesInH(g))
        all_zero <- AllZero(g)
        i <- i + 1
    }
    
    return(global_l)
}

GraphInit <- function(g) {

    V(g)$id <- 1:length(V(g))
#    g <- LinkageFunctionInit(g)
    return(g)
}

# algorithm to find the best cluster from the current node selection
MultipartiteBestCluster <- function(g) {

    g <- LinkageFunctionInit(g)
    Gamma <- c()
    FGamma <- -INFTY
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
    
    l <- list("solution" = Gamma, "value" = FGamma)
    return(list("graph" = g, "cluster" = l)) 
}

AllZero <- function(g) {

    zero_val <- which(LinkageVal(g, GetNodesInH(g)) == 0)
    return(length(zero_val) == length(GetNodesInH(g)))
}

# return a vector where the components are the values of the linkage function computed on each node of the graph
LinkageFunctionInit <- function(g) {

    # initialise the optimal set
    V(g)$is.inH <- rep(TRUE, length(V(g)))
    # initialise the linkage value = sum of the weights of incident edges, for every node
    V(g)$linkage.value <- sapply(V(g), function(x) sum(edge_attr(g, "weight", incident(g,x))))
    
    return(g)
}

# remove the node indexes (1-based) and update the linkage function
LinkageFunctionUpdate <- function(g, node.id) {

    node.idx <- Idx2Id(g, node.id)
    V(g)$is.inH[node.idx] <- FALSE
    edges1 <- unlist(incident_edges(g, node.idx)) # edges incident in node.idx
    edges2 <- unlist(incident_edges(g, V(g)[V(g)$is.inH])) # edges incident in any node of H
    edges <- edges1[edges1 %in% edges2] # edges connecting node.idx and all the nodes in H
    nodes <- c(ends(g, edges))
    nodes <- nodes[!(nodes %in% node.idx)] # nodes adjacent to node.idx, in the same order as the edges, with repeats
    val_edges <- edge_attr(g, "weight", edges) # weight of the edges incident to the nodes above
    idx <- unique(nodes)
    val <- sapply(idx, function(x) sum(val_edges[nodes == x])) # each entry is the sum for the nodes in H
    V(g)$linkage.value[idx] <- V(g)$linkage.value[idx] - 2*val

    return(g)
}

# find the node with the minimum linkage value
FindArgMinLinkage <- function(g) {

    val <- min(V(g)$linkage.value[V(g)$is.inH])
    idx <- as_ids(V(g)[V(g)$linkage.value <= val & V(g)$is.inH]) 
    return(Idx2Id(g, idx))
}

LinkageVal <- function(g, node.id) {

    node.idx <- Id2Idx(g, node.id)
    v <- c(node.idx)
    return(V(g)$linkage.value[node.idx])
}

DeleteNodes <- function(g, node.id) {

    node.idx <- Id2Idx(g, node.id)
    g <- delete_vertices(g, node.idx)
    return(g)
}

GetNodesInH <- function(g) {

    idx <- as_ids(V(g)[V(g)$is.inH])
    return(Idx2Id(g, idx))
}

SetNodesInH <- function(g, node.id) {

    idx <- Id2Idx(g, node.id)
    V(g)$is.inH <- as_ids(V(g)) %in% idx
    return(g)
}

Id2Idx <- function(g, node.id) {

    if (sum(V(g)$id %in% node.id) < length(node.id)) {
        stop(paste("Node id", node.id, "not found in this graph\n"))
    }
    return(as_ids(V(g)[match(node.id,V(g)$id)]))
}

Idx2Id <- function(g, node.idx) {

    if (length(node.idx) == 0) {
        return(NULL)
    }
    if (max(node.idx) > length(V(g))) {
        stop(paste("Node idx", node.idx, "exceeds graph size\n"))
    }
    return(V(g)$id[node.idx])
}


