###### Convert matrix list to big matrix ######
mat.convert <- function(mm){
    m <- dim(mm)[1]
    res <- matrix(list(), nrow = m, ncol = 1)
    for(i in 1:m){
        res[[i,1]] <- mm[[i,1]]
    }
    for(i in 1:m){
        for(j in 2:m){
            res[[i,1]] <- cbind(res[[i,1]], mm[[i,j]])
        }
    }
    result <- res[[1,1]]
    for(i in 2:m){
        result <- rbind(result, res[[i,1]])
    }
    return(as.matrix(result))
}
