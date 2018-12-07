###### ADMM solver ######
mod.admm <- function(X, Y, D, S, Lambda, lambda, rho_lambda, rho, tol){
    ### S is location structure pattern matrix
    ### Lambda is precision matrix list
    ### alpha is relaxation parameter
    ### lambda is tuning parameter vector
    ### rho is converge stepsize
    e <- dim(D)[1]
    p <- dim(D)[2]
    len <- dim(X[[1]])[1] ### Used time series length
    Gamma = Zeta <- matrix(runif(m^2), m, m)
    Zeta.old = Eta <- matrix(0, m, m)
    Delta = B = C = Omega1 = Omega2 <- matrix(list(), m, m)
    for(i in 1:m){
        for(j in 1:m){
            Delta[[i,j]] = B[[i,j]] = matrix(runif(p^2), p, p)
            Omega1[[i,j]] <- matrix(0, p, p)
            C[[i,j]] <- D %*% Delta[[i,j]]
            Omega2[[i,j]] <- matrix(0, e, p)
        }
    }
    B.old = C.old <- matrix(list(), m, m)
    dDelta <- matrix(list(), m, m)
    
    ### Recording Primal and Dual residuals
    Prim.gamma = Dual.gamma <- matrix(0, m, m)
    Prim.delta = Dual.delta <- matrix(list(), m, m)
    Prim.fusion = Dual.fusion <- matrix(list(), m, m)
    for(i in 1:m){
        for(j in 1:m){
            Prim.delta[[i,j]] = Dual.delta[[i,j]] <- matrix(0, p, p)
            Prim.fusion[[i,j]] = Dual.fusion[[i,j]] <- matrix(0, p, p)
        }
    }
    
    obj.result <- c()
    Prim.res <- 1
    Dual.res <- 1
    
    obj.update <- 0
    obj <- objective(X, Y, Gamma, Delta, Lambda, Zeta, B, C, lambda, rho_lambda)
    
    tol.pri <- sqrt(p) * tol
    tol.dual <- sqrt(len) * tol
    rep <- 1
    while((Prim.res >= tol.pri) || (Dual.res >= tol.dual)){
        ### Update Gamma by elements
        for(i in 1:m){
            for(j in 1:m){
                ### Recording some essential values
                y <- correspond.sum(X, i, j, Zeta, B)
                q <- (tr(Lambda[[i]] %*% t(X[[j]] %*% t(B[[i,j]])) %*% (Y[[i]] - y) / len) + rho * (Zeta[i,j] - Eta[i,j]))
                
                ### Update Gamma[i,j]
                Gamma[i,j] <- q / (tr(Lambda[[i]] %*% t(X[[j]] %*% t(B[[i,j]])) %*% (X[[j]] %*% t(B[[i,j]])) / len) + rho)
                
                ### Update Zeta[i,j]
                Zeta.old[i,j] <- Zeta[i,j]
                Zeta[i,j] <- soft.threshold.scalar(Gamma[i,j] + Eta[i,j], lambda[1] / rho)
                
                ### Update Eta[i,j]
                Eta[i,j] <- Eta[i,j] + (Gamma[i,j] - Zeta[i,j])
            }
        }
        ### Primal and Dual residuals
        Prim.gamma <- Gamma - Zeta
        Dual.gamma <- rho * (Zeta - Zeta.old)
        
        ### Update Delta by blocks
        delta.norm = ddelta.norm <- 0
        B.norm = C.norm <- 0
        omega1.norm = omega2.norm <- 0
        for(i in 1:m){
            for(j in 1:m){
                U <- factorization(X[[j]], D, Lambda[[i]], rho, Zeta[i,j], len)$U
                L <- factorization(X[[j]], D, Lambda[[i]], rho, Zeta[i,j], len)$L
                
                ### Only consider skinny matrix, recording essential values
                y <- correspond.sum(X, i, j, Zeta, B)
                q <- (Zeta[i,j] * Lambda[[i]] %*% t(X[[j]]) %*% (Y[[i]] - y) / len) + (rho * (B[[i,j]] - Omega1[[i,j]]) + rho * t(D) %*% (C[[i,j]] - Omega2[[i,j]]))
                
                ### Update Delta
                Delta[[i,j]] <- solve(U, solve(L, q))
                
                ### Update B
                B.old[[i,j]] <- B[[i,j]] 
                B[[i,j]] <- soft.threshold(Delta[[i,j]] + Omega1[[i,j]], lambda[2] / rho)
                
                ### Update C
                C.old[[i,j]] <- C[[i,j]]
                dDelta[[i,j]] <- D %*% Delta[[i,j]]
                C[[i,j]] <- soft.threshold(dDelta[[i,j]] + Omega2[[i,j]], lambda[3] / rho)
                
                ### Update Omega's
                Omega1[[i,j]] <- Omega1[[i,j]] + (Delta[[i,j]] - B[[i,j]])
                Omega2[[i,j]] <- Omega2[[i,j]] + (dDelta[[i,j]] - C[[i,j]])
                
                ### Check the spatial structure pattern with matrix S
                B[[i,j]][(B[[i,j]] != 0 & S == 0)] <- 0
                
                ### Recording norms for updating stop criterion
                delta.norm <- max(delta.norm, norm(Delta[[i,j]], "F"))
                ddelta.norm <- max(ddelta.norm, norm(D %*% Delta[[i,j]], "F"))
                B.norm <- max(B.norm, norm(B[[i,j]], "F"))
                C.norm <- max(C.norm, norm(C[[i,j]], "F"))
                omega1.norm <- (omega1.norm + norm(Omega1[[i,j]], "F"))
                omega2.norm <- (omega2.norm + norm(t(D) %*% Omega2[[i,j]], "F"))
                
                ### Recording primal and dual residuals for Delta and Fusion
                Prim.delta[[i,j]] <- Delta[[i,j]] - B[[i,j]]
                Prim.fusion[[i,j]] <- D %*% Delta[[i,j]] - C[[i,j]]
                Dual.delta[[i,j]] <- rho * (B[[i,j]] - B.old[[i,j]])
                Dual.fusion[[i,j]] <- rho * t(D) %*% (C[[i,j]] - C.old[[i,j]])
            }
        }
        Prim.delta.mat <- mat.convert(Prim.delta)
        Prim.fusion.mat <- mat.convert(Prim.fusion)
        Dual.delta.mat <- mat.convert(Dual.delta)
        Dual.fusion.mat <- mat.convert(Dual.fusion)
        
        Prim.res <- norm(Prim.gamma, "F") + norm(Prim.delta.mat, "F") + norm(Prim.fusion.mat, "F")
        Dual.res <- norm(Dual.gamma, "F") + norm(Dual.delta.mat, "F") + norm(Dual.fusion.mat, "F")
        
        ### Update stepsize by residuals
        if(Prim.res > 10 * Dual.res){
            rho <- 2 * rho
        }else if(Prim.res < 0.1 * Dual.res){
            rho <- rho / 2
        }else{
            rho <- rho
        }
        obj.update <- objective(X, Y, Gamma, Delta, Lambda, Zeta, B, C, lambda, rho_lambda)
        tol.pri <- sqrt(p) * tol + 1e-3 * max(norm(Gamma, "F"), norm(Zeta, "F"), delta.norm, ddelta.norm, B.norm, C.norm)
        tol.dual <- sqrt(len) * tol + 1e-3 * max(norm(Eta, "F"), omega1.norm, omega2.norm)
        obj.result <- c(obj.result, obj.update)
        #print(B[[i,j]][1:5,1:5])
        print(c(rep, obj.update, rho, Prim.res, Dual.res, tol.pri, tol.dual))
        rep <- rep + 1
    }
    A.est <- matrix(list(), m, m)
    for(i in 1:m){
        for(j in 1:m){
            A.est[[i,j]] <- Zeta[i,j] * B[[i,j]]
        }
    }
    return(list(Gamma.est = Zeta, Delta.est = B, A.est = A.est, obj.values = obj.update, obj.seq = obj.result))
}
