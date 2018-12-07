###### Sum function ######
correspond.sum <- function(X, i, j, Gamma, Delta){
    y <- 0
    for(k in 1:m){
        y <- y + Gamma[i,k] * X[[k]] %*% t(Delta[[i,k]])
    }
    return(y - Gamma[i,j] * X[[j]] %*% t(Delta[[i,j]]))
}
