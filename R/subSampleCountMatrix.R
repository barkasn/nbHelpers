


#' Subsample a count matrix to a specific number of reads
#' @param mat sparse matrix of class dgCMatrix
#' @param final.count final number of reads the matrix should have
#' @param verbose logical verbosity
#' @return a new dgCMatrix subsample to the requested number of reads
#' @export subSampleCountMatrix
subSampleCountMatrix <-
  function (mat,
            final.count = 1e+07,
            verbose = FALSE)
  {
    if (class(mat) != "dgCMatrix") {
      stop("mat is not of class dgCMatrix")
    }
    final.count <- floor(final.count)
    i <- 1

    while (sum(mat@x) > final.count) {
      if (verbose) {
        cat(paste0('Iteration ', i, '... mat sum ', sum(mat@x), '\n'))

      }
      i <- i + 1

      mat <- drop0(mat)
      current.count <- sum(mat@x)
      if (current.count < final.count) {
        stop("Reads in matrix are fewer than the final counts")
      }
      reads.remove <- current.count - final.count
      pos.rem <- sample(length(mat@x), reads.remove, replace = TRUE)
      pos.rem.tbl <- table(pos.rem)
      pos.rem.num <- as.numeric(pos.rem.tbl)
      names(pos.rem.num) <- names(pos.rem.tbl)
      x <- mat@x
      rmvec <- pos.rem.num[1:length(x)]
      rmvec[is.na(rmvec)] <- 0
      x <- x - rmvec
      x <- pmax(x, 0)
      mat@x <- x
      mat <- Matrix::drop0(mat)
    }
    mat
  }


#' Subsample a count matrix to a specific number of reads, improved version
#' ensures that every read has the same probability of being dropped
#' @param mat sparse matrix of class dgCMatrix
#' @param final.count final number of reads the matrix should have
#' @param seed seed for random number generator
#' @return a new dgCMatrix subsample to the requested number of reads
#' @export subSampleCountMatrix2
subSampleCountMatrix2 <- function (mat,final.count = 1e+06, seed=NULL)
{
  if (class(mat) != "dgCMatrix") {
    stop("mat is not of class dgCMatrix")
  }
  ##
  if (!is.null(seed)) {
    set.seed(seed)
  }
  ##
  final.count <- floor(final.count)
  ## remove empty positions
  mat <- drop0(mat)
  ## calculate how many reads to remove
  current.count <- sum(mat@x)
  if (current.count < final.count) {
    stop("Reads in matrix are fewer than the final counts")
  }
  reads.remove <- current.count - final.count
  ## Sample from the total number of reads
  srr <- sample(sum(mat@x),reads.remove, replace=FALSE)
  ## map the read ids to positions on the x vector
  csx <- cumsum(mat@x)
  reads.remove.entries <- table(findInterval(srr, csx))
  p <- reads.remove.entries[1:length(mat@x)]
  p[is.na(p)] <- 0
  p <- as.numeric(p)
  ## remove the read
  mat@x <- mat@x - p
  ## drop any newly appearing 0s
  mat <- Matrix::drop0(mat)
}

