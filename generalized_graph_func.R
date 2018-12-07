###### Selecting bridges from graph communities ######
select.bridge <- function(g, eb){
    D <- as.matrix(getDgSparse(g))
    n.communities <- length(communities(eb))
    n.edges <- dim(D)[1]
    
    coor <- matrix(0, nrow = n.edges, 2)
    for(i in 1:n.edges){
        coor[i,] <- which(D[i,] != 0)
    }
    bridge <- NULL
    for(i in 1:n.edges){
        for(j in 1:n.communities){
            if((coor[i,1] %in% eb[[j]]) && !(coor[i,2] %in% eb[[j]])){
                for(k in 1:n.communities){
                    if(coor[i,2] %in% eb[[k]]){
                        bridge <- rbind(bridge, c(coor[i,], i, j, k))
                    }
                }
            }
        }
    }
    return(list(bridge = bridge, D = D))
}

###### Weighting function ######
weight.mat <- function(bridges.result, eb){
    bridges <- bridges.result$bridge
    D <- bridges.result$D
    n.bridge <- dim(bridges)[1]
    n.edge <- dim(D)[1]
    n.communities <- length(communities(eb))
    
    n.vex_in_com <- c()
    for(i in 1:n.communities){
        n.vex_in_com[i] <- length(communities(eb)[[i]])
    }
    
    ### Compute bridge weights
    if(!is.null(bridges)){
        for(i in 1:n.bridge){
            row <- bridges[i,3]
            col.from <- bridges[i,1]
            col.to <- bridges[i,2]
            community.from <- bridges[i,4]
            community.to <- bridges[i,5]
            
            weight <- n.vex_in_com[community.to] / (n.vex_in_com[community.from] + n.vex_in_com[community.to])
            D[row, col.from] <- -weight
            D[row, col.to] <- weight
        }
    }
    
    ### Compute local edge weights
    for(i in 1:n.edge){
        v.from <- which(D[i,] == -1)
        v.to <- which(D[i,] == 1)
        if((length(v.from) != 0) && (length(v.to) != 0)){
            for(j in 1:n.communities){
                if((v.from %in% eb[[j]]) && (v.to %in% eb[[j]])){
                    neighbor.from <- neighborhood(g, 1, v.from)[[1]]
                    neighbor.to <- neighborhood(g, 1, v.to)[[1]]
                    
                    weight <- length(intersect(neighbor.to, neighbor.from)) / length(union(neighbor.from, neighbor.to))
                    #print(c(v.from, v.to, weight))
                    D[i,v.from] <- -weight
                    D[i,v.to] <- weight
                }
            }
        }
    }
    return(D)
}

###### Generalized spatial structure functions ######
edges.select <- function(dmat, edge.coef = 0.17, bridge.coef = 0.22){
    p <- dim(dmat)[1]
    dist.max <- max(apply(dmat, 2, max, na.rm = TRUE))
    dist.edge <- edge.coef * dist.max
    dist.bridge <- bridge.coef * dist.max
    edges <- c()
    bridges <- c()
    for(i in 1:(p-1)){
        for(j in (i+1):p){
            if(dmat[i,j] <= dist.edge){
                edges <- c(edges, dmat[i,j])
            }else if(dmat[i,j] > dist.edge && (dmat[i,j] <= dist.bridge)){
                bridges <- c(bridges, dmat[i,j])
            }
        }
    }
    return(list(edges = edges, bridges = bridges))
}

weighted.incident.mat <- function(dmat, edges.result){
    edges <- edges.result$edges
    bridges <- edges.result$bridges
    D.edges <- matrix(0, nrow = length(edges), ncol = p)
    D.bridges <- matrix(0, nrow = length(bridges), ncol = p)
    for(i in 1:length(edges)){
        id.edge <- which(dmat == edges[i], arr.ind = TRUE)
        pt.from <- id.edge[1,2]
        pt.to <- id.edge[1,1]
        D.edges[i, pt.from] <- -1 / edges[i]
        D.edges[i, pt.to] <- 1 / edges[i]
    }
    for(i in 1:length(bridges)){
        id.bridge <- which(dmat == bridges[i], arr.ind = TRUE)
        pt.from <- id.bridge[1,2]
        pt.to <- id.bridge[1,1]
        D.bridges[i, pt.from] <- -1 / bridges[i]
        D.bridges[i, pt.to] <- 1 / bridges[i]
    }
    D <- rbind(D.edges, D.bridges)
    return(D)
}

spatial.structure <- function(dmat, edges.result){
    edges <- edges.result$edges
    bridges <- edges.result$bridges
    adjacency.mat <- matrix(0, p, p)
    for(i in 1:length(edges)){
        index <- which(dmat == edges[i], arr.ind = TRUE)
        adjacency.mat[index[1,1], index[1,2]] <- 1
        adjacency.mat[index[2,1], index[2,2]] <- 1
    }
    for(i in 1:length(bridges)){
        index <- which(dmat == bridges[i], arr.ind = TRUE)
        adjacency.mat[index[1,1], index[1,2]] <- 1
        adjacency.mat[index[2,1], index[2,2]] <- 1
    }
    spatial.mat <- adjacency.mat + diag(p)
    return(list(adjacency.mat = adjacency.mat, spatial.mat = spatial.mat))
}

generate.graph <- function(spatial.result){
    adjacency.mat <- spatial.result$adjacency.mat
    g <- graph_from_adjacency_matrix(adjacency.mat, "undirected")
    eb <- cluster_edge_betweenness(g)
    return(list(g = g, eb = eb))
}

plot.coor <- function(coor.mat, dmat, edges.result){
    edges <- edges.result$edges
    bridges <- edges.result$bridges
    plot(coor.mat)
    for(i in 1:length(edges)){
        id.edge <- which(dmat == edges[i], arr.ind = TRUE)
        lines(coor.mat[id.edge[1,],], col = "blue")
    }
    for(i in 1:length(bridges)){
        id.bridge <- which(dmat == bridges[i], arr.ind = TRUE)
        lines(coor.mat[id.bridge[1,],], col = "red")
    }
}
