###### Objective function ######
objective <- function(X, Y, Gamma, Delta, Lambda, Zeta, B, C, lambda, rho_lambda){
    m <- dim(Zeta)[1]
    p <- dim(B[[1,1]])[1]
    len <- dim(X[[1]])[1]
    lik <- 0
    pen <- 0
    ### Define a lambda matrix to modifty incident matrix
    for(i in 1:m){
        partial.sum <- matrix(0, nrow = len, ncol = p)
        for(j in 1:m){
            partial.sum <- partial.sum + Gamma[i,j] * X[[i]] %*% t(Delta[[i,j]])
        }
        lik <- lik + tr(Lambda[[i]] %*% t(Y[[i]] - partial.sum) %*% (Y[[i]] - partial.sum)) / len
    }
    pen.delta = pen.precision <- 0
    for(i in 1:m){
        for(j in 1:m){
            pen.delta <- pen.delta + lambda[2] * norm(B[[i,j]], "o") + lambda[3] * norm(C[[i,j]], "o")
        }
        pen.precision <- pen.precision + rho_lambda * norm(Lambda[[i]], "o")
    }
    pen <- pen.delta + lambda[1] * norm(Zeta, "o") + pen.precision
    return(0.5 * lik + pen)
}
