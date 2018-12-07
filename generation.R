###### Generating data function ######
generate.data <- function(p, m, n, transition.matrix, variance.mat, burnin = 50){
    data <- matrix(0, nrow = n, ncol = m*p)
    data[1,] <- rnorm(m*p)
    for(t in 2:n){
        error <- c()
        for(i in 1:m){
            err <- as.matrix(mvrnorm(1, rep(0, p), variance.mat[[i,1]]))
            error <- rbind(err, error)
        }
        data[t,] <- data[t-1,] %*% t(transition.matrix) + t(error)
    }
    X = Y <- list()
    for(i in 1:m){
        X[[i]] <- data[burnin:(n-1), ((i-1)*p+1):(i*p)]
        Y[[i]] <- data[(burnin + 1):n, ((i-1)*p+1):(i*p)]
    }
    return(list(X = X, Y = Y))
}
