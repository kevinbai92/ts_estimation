###### BIC criterion function ######
estimated.par <- function(fit){
    ### since we have only 3 tuning parameters. It is a 3-loop.
    flag.gamma = flag.A <- 0
    for(i in 1:m){
        for(j in 1:m){
            if(fit$Gamma.est[i,j] != 0){
                flag.gamma <- flag.gamma + 1
            }
            for(k in 1:p){
                for(l in 1:p){
                    if(fit$A.est[[i,j]][k,l] != 0){
                        flag.A <- flag.A + 1
                    }
                }
            }
        }
    }
    return(flag.A + flag.gamma)
}

bic <- function(n, k, obj.value){
    return(log(n)*k - 2*log(obj.value))
}
