
#' mix cells in a matrix with a common background profile generated from all the cell
#' @param matrix of genes x cells
#' @param bg.prog proportion of read to come from the background
#' @param n.cores number of cores t use
#' @param verbose print out dots to show while processing
#' @return object of class dgCMatrix
#' @export mixCellsInMatrix
mixCellsInMatrix <- function(mat, bg.prop, n.cores=1, verbose=FALSE) {
  ## calc the average distribution of reads (the background)
  avg.dist <- apply(mat,1,sum)
  ## for every cell
  new.mat <- do.call(cbind, mclapply(1:ncol(mat), function(ci) {
    ## pick these many reads from one of the other distributions
    cell.expr.v <- mat[,ci]
    reads.in.cell <- sum(cell.expr.v)
    ## calculate reads from each pool
    reads.bg.n <- floor(reads.in.cell * bg.prop)
    reads.cell.n <- ceiling(reads.in.cell * (1 - bg.prop))
    ## sample the background
    srr.bg <- sample(sum(avg.dist),reads.bg.n, replace=FALSE)
    csx.bg <- cumsum(avg.dist)
    toi.bg <- table(findInterval(srr.bg,csx.bg,rightmost.closed = FALSE, left.open = TRUE) + 1)
    x <- as.numeric(toi.bg); names(x) <- names(toi.bg); toi.bg <-x; rm(x)
    toi.bg.2 <- toi.bg[match(as.character(1:length(avg.dist)), names(toi.bg))]
    toi.bg.2[is.na(toi.bg.2)] <- 0
    ## sample the cell
    srr.cl <- sample(sum(cell.expr.v), reads.cell.n, replace=FALSE)
    csx.cl <- cumsum(cell.expr.v)
    toi.cl <- table(findInterval(srr.cl, csx.cl, rightmost.closed = FALSE, left.open = TRUE) + 1)
    x <- as.numeric(toi.cl); names(x) <- names(toi.cl); toi.cl <- x; rm(x)
    toi.cl.2 <- toi.cl[match(as.character(1:length(cell.expr.v)), names(toi.cl))]
    toi.cl.2[is.na(toi.cl.2)] <- 0
    ## put together
    final.cv <- toi.bg.2 + toi.cl.2
    names(final.cv) <- NULL
    ## Sanity check
    ##sum(abs(final.cv - cell.expr.v))
    ## return
    if(verbose) { if (ci %% 100 == 0) { cat('.') } }
    final.cv
  }, mc.cores= n.cores))
  ## copy names
  colnames(new.mat) <- colnames(mat)
  rownames(new.mat) <- rownames(mat)
  ## return  as a sparse matrix
  as(new.mat, 'dgCMatrix')
}
