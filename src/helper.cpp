#include "helper.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat matrixCorr(arma::mat X, arma::mat Y) {
    arma::mat c=arma::cor(X,Y);
    return c;
}
