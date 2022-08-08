#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_Psi_hat(int N, List parameters) {
  NumericMatrix Psi(N, N*N*N);
  NumericVector gamma = parameters["gamma"];
  NumericVector scatter = parameters["scatter"];
  List a = parameters["a"];
  double a41 = Rcpp::as<Rcpp::List>(parameters["a"])["a41"];
  double a42 = Rcpp::as<Rcpp::List>(parameters["a"])["a42"];
  double a43 = Rcpp::as<Rcpp::List>(parameters["a"])["a43"];
  for(int i=0; i < N;i++) {
    for(int j=0; j < N;j++) {
      for(int k=0; k < N;k++){
        for(int l=0; l < N;l++) {
          Psi(i, j* N*N + k*N + l) = a41 * gamma(i) * gamma(j) * gamma(k) * gamma(l) +
            (a42/6) * (gamma(i) * gamma(j) * scatter(k,l) + gamma(i) * gamma(k) * scatter(j,l) + gamma(i) * gamma(l) * scatter(j,k) + gamma(j) * gamma(k) * scatter(i,l) + gamma(j) * gamma(l) * scatter(i,k) + gamma(k) * gamma(l) * scatter(i,j) ) +
            (a43/3) * (scatter(i,j) * scatter(k,l) + scatter(i,k) * scatter(j,l) + scatter(i,l) * scatter(j,k));
        }
      }
    }
  }
  return(Psi);
}
 
