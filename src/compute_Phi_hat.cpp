#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_Phi_hat(int N, List parameters) {
  NumericMatrix Phi(N, N*N);
  NumericVector gamma = parameters["gamma"];
  NumericVector scatter = parameters["scatter"];
  List a = parameters["a"];
  double a31 = Rcpp::as<Rcpp::List>(parameters["a"])["a31"];
  double a32 = Rcpp::as<Rcpp::List>(parameters["a"])["a32"];
  for(int i=0; i < N;i++) {
    for(int j=0; j < N;j++) {
      for(int k=0; k < N;k++){
        Phi(i, j* N + k) = a31 * gamma(i) * gamma(j) * gamma(k) +
          (a32/3) * (gamma(i) * scatter(j,k) + gamma(j) * scatter(i,k) + gamma(k) * scatter(i,j) );
      }
    }
  }
  return(Phi);
}

