// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_Phi_hat
NumericMatrix compute_Phi_hat(int N, List parameters);
RcppExport SEXP _ParametricHighOrderPortfolios_compute_Phi_hat(SEXP NSEXP, SEXP parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type parameters(parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Phi_hat(N, parameters));
    return rcpp_result_gen;
END_RCPP
}
// compute_Psi_hat
NumericMatrix compute_Psi_hat(int N, List parameters);
RcppExport SEXP _ParametricHighOrderPortfolios_compute_Psi_hat(SEXP NSEXP, SEXP parametersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type parameters(parametersSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Psi_hat(N, parameters));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ParametricHighOrderPortfolios_compute_Phi_hat", (DL_FUNC) &_ParametricHighOrderPortfolios_compute_Phi_hat, 2},
    {"_ParametricHighOrderPortfolios_compute_Psi_hat", (DL_FUNC) &_ParametricHighOrderPortfolios_compute_Psi_hat, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ParametricHighOrderPortfolios(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
