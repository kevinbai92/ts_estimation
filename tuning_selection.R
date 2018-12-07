###### Choosing tuning range ######
lambda.range <- function(from, to, len, S, Lambda, rho_lambda, rho, tol){
    range <- array(0, dim = c(rep(len, 3)))
    range.seq.gamma <- sort(runif(len+1, from[1], to[1]))
    range.seq.delta <- sort(runif(len+1, from[2], to[2]))
    range.seq.fusion <- sort(runif(len+1, from[3], to[3]))
    for(i in 1:len){
        for(j in 1:len){
            for(k in 1:len){
                lambda.gamma <- range.seq.gamma[i]
                lambda.delta <- range.seq.delta[j]
                lambda.fusion <- range.seq.fusion[k]
                
                fit <- mod.admm(X, Y, D, S, Lambda, lambda = c(lambda.gamma, lambda.delta, lambda.fusion), rho_lambda, rho = rho, tol = tol)
                trans.mat <- mat.convert(fit$A.est)
                npar <- estimated.par(fit)
                range[i,j,k] <- bic(n, npar, fit$obj.values)
            }
        }
    }
    coor.est <- which(range == min(range, na.rm = TRUE), arr.ind = TRUE)
    new.from <- c(range.seq.gamma[coor.est[1]], range.seq.delta[coor.est[2]], range.seq.fusion[coor.est[3]])
    new.to <- c(range.seq.gamma[coor.est[1]+1], range.seq.delta[coor.est[2]+1], range.seq.fusion[coor.est[3]+1])
    return(list(from = new.from, to = new.to))
}

###### Choosing tuning parameters ######
lambda.select <- function(from, to, lambda.len, S, Lambda, rho_lambda, rho, tol){
    ### 'from' and 'to' are both range vectors for tuning parameters gamma, delta and fusion
    bic.result <- array(0, dim = c(rep(lambda.len, 3)))
    lambda.seq.gamma <- runif(lambda.len, from[1], to[1])
    lambda.seq.delta <- runif(lambda.len, from[2], to[2])
    lambda.seq.fusion <- runif(lambda.len, from[3], to[3])
    
    for(i in 1:lambda.len){
        for(j in 1:lambda.len){
            for(k in 1:lambda.len){
                lambda.gamma <- lambda.seq.gamma[i]
                lambda.delta <- lambda.seq.delta[j]
                lambda.fusion <- lambda.seq.fusion[k]
                
                fit <- mod.admm(X, Y, D, S, Lambda, lambda = c(lambda.gamma, lambda.delta, lambda.fusion), rho_lambda, rho = rho, tol = tol)
                npar <- estimated.par(fit)
                bic.result[i,j,k] <- bic(n, npar, fit$obj.values)
            }
        }
    }
    lambda.est <- c(lambda.seq.gamma[which(bic.result == min(bic.result, na.rm = TRUE), arr.ind = TRUE)[1]], 
                    lambda.seq.delta[which(bic.result == min(bic.result, na.rm = TRUE), arr.ind = TRUE)[2]],
                    lambda.seq.fusion[which(bic.result == min(bic.result, na.rm = TRUE), arr.ind = TRUE)[3]])
    return(lambda.est)
}
