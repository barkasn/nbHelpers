#' @importFrom utils read.table
#' @importFrom Matrix readMM
NULL

#' Generates a random positive definite matrix
#' @description generates a random positive definite matrix with the specified
#' eigenvalues. Code from: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.html
#' @param n dimensions
#' @param ev optional eigenvalues
#' @return a random positive definite matrix with requested eigenvalues
#' @export randomPositiveDefiniteMatrix
randomPositiveDefiniteMatrix <- function(n, ev = runif(n, 0, 10)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

#' Generates data with covariance structure specified by the given covariance matrix
#' @description Generates the requested number of observations that follow the given
#' covariance structure specified by the given matrix
#' @param covM a posititive definite covariance matrix
#' @param nobs the number of observations requestest
#' @return data.frame object with the observations as rows and the variables as columns
#' @import matrixcalc
#' @export getDataWithCovarianceMatrix
getDataWithCovarianceMatrix <- function(covM, nobs) {
  # https://www.r-bloggers.com/simulating-data-following-a-given-covariance-structure/

  covDim <- dim(covM)
  nvars <- covDim[2]

  if (nobs < 1)
    stop('No observations were requested')

  if (covDim[1] != covDim[2])
    stop('Covariance matrix is not square!')

  if (!matrixcalc::is.positive.definite(covM)) {
    stop('Covariance matrix is not positive definite!')
  }

  as.data.frame(t(t(chol(covM)) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)))
}
