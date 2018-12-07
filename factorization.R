###### Factorization function ######
factorization <- function(A, D, Lambda, rho, gamma, n){
    p <- dim(A)[2]
    ### Only for skinny matrix 
    U <- chol(gamma^2 * Lambda %*% t(A) %*% A / n + rho * (t(D) %*% D + diag(p)))
    L <- t(U)
    return(list(U = U, L = L))
}