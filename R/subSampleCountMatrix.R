
#' Subsample a count matrix to a specific number of reads
#' @param mat sparse matrix of class dgCMatrix
#' @param final.count final number of reads the matrix should have
#' @return a new dgCMatrix subsample to the requested number of reads
#' @export subSampleCountMatrix
subSampleCountMatrix <- function(mat, final.count=1e7) {
    if(class(mat) != 'dgCMatrix') {
        stop('mat is not of class dgCMatrix')
    }
    mat <- drop0(mat)
    ## Get current count and check
    current.count <- sum(mat@x)
    if (current.count < final.count) {
        stop('Reads in matrix are fewer than the final counts');
    }
    ## reads to remove
    reads.remove <- current.count - final.count
    ## sample positions to remove
    pos.rem <- sample(length(mat@x),reads.remove,replace=TRUE)
    pos.rem.tbl <- table(pos.rem)
    ## convert to numeric vector
    pos.rem.num <- as.numeric(pos.rem.tbl)
    names(pos.rem.num) <- names(pos.rem.tbl)
    ## get the x vector
    x <- mat@x
    ## construct remove vector
    rmvec <- pos.rem.num[1:length(x)]
    rmvec[is.na(rmvec)] <- 0
    ## remove reads 
    x <- x - rmvec
    x <- pmax(x, 0)
    ## put matrix back together
    mat@x <- x
    mat <- Matrix::drop0(mat)
    ## due to negatives we have probably not hit the target, iterate
    if(sum(mat@x) > final.count) {
        mat <- subSampleCountMatrix(mat,final.count)
    }
    ## return
    mat
}

