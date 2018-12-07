###### Trace function ######
tr <- function(M){
    if(ncol(M) != nrow(M)){
        print("Not a square matrix")
    }else{
        return(sum(diag(M)))
    }
}
