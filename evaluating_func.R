###### Evaluation methods ######
classification <- function(M.est, M.true){
    m <- dim(M.est)[1]
    n <- dim(M.est)[2]
    tp = tn <- 0
    fp = fn <- 0
    for(i in 1:m){
        for(j in 1:n){
            if(M.true[i,j] != 0){
                if(M.est[i,j] != 0){
                    tp <- tp + 1
                }else if(M.est[i,j] == 0){
                    fn <- fn + 1
                }
            }else if(M.true[i,j] == 0){
                if(M.est[i,j] != 0){
                    fp <- fp + 1
                }else if(M.est[i,j] == 0){
                    tn <- tn + 1
                }
            }
        }
    }
    return(list(tp = tp, fp = fp, tn = tn, fn = fn))
}

###### Evaluaton functions ######
eval.func <- function(classify, method){
    tp <- classify$tp
    fp <- classify$fp
    tn <- classify$tn
    fn <- classify$fn
    if(method == "SEN"){
        return(tp / (tp + fn))
    }else if(method == "SPE"){
        return(tn / (tn + fp))
    }else if(method == "PRE"){
        return(tp / (tp + fp))
    }else if(method == "Fscore"){
        f <- 2*tp / (2*tp + fp + fn)
        return(f)
    }else if(method == "ACC"){
        acc <- (tp + tn) / (tp + tn + fp + fn)
        return(acc)
    }else if(method == "MCC"){
        mcc <- (tp * tn - fp * fn) / sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
        return(mcc)
    }
}

re.eval <- function(M.est, M.true){
    if(norm(M.true, "F") != 0){
        return(norm(M.est - M.true, "F") / norm(M.true, "F"))
    }else{
        return(0)
    }
}
