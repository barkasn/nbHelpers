// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// matrixCorr
arma::mat matrixCorr(arma::mat X, arma::mat Y);
RcppExport SEXP _nbHelpers_matrixCorr(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixCorr(X, Y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nbHelpers_matrixCorr", (DL_FUNC) &_nbHelpers_matrixCorr, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_nbHelpers(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}