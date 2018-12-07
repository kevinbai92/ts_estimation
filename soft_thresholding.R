###### Soft-thresholding functions ######
soft.threshold <- function(A, lambda){
    m <- dim(A)[1]
    n <- dim(A)[2]
    Z <- matrix(0, nrow = m, ncol = n)
    for(i in 1:m){
        for(j in 1:n){
            Z[i,j] <- (max(0, A[i,j] - lambda) - max(0, -A[i,j] - lambda))
        }
    }
    return(Z)
}

soft.threshold.scalar <- function(a, lambda){
    r <- max(0, a - lambda) - max(0, -a - lambda)
    if(r < 0){
        return(0)
    }else{
        return(r)
    }
}
