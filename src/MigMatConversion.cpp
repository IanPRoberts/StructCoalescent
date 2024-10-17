// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>

using namespace Rcpp;

//' @title Forward-in-time migration matrix conversion
//' @description Computes the forward-in-time migration matrix associated with a given backwards-in-time migration matrix, coalescent rate and time scale
//' @param bit_mm NumericMatrix Backward-in-time migration matrix
//' @param coal_rate NumericVector Coalescent rates
//' @returns NumericMatrix Forward-in-time migration matrix
//'
//' @export
//'
// [[Rcpp::export]]

NumericMatrix FitMigMatC(NumericMatrix bit_mm, NumericVector coal_rate) {
  int n_deme = bit_mm.nrow();
  NumericMatrix out(n_deme, n_deme);

  for (int i = 0; i < n_deme; ++i){
    for (int j = 0; j < n_deme; ++j){
      out(i,j) = bit_mm(j,i) * coal_rate(i) / coal_rate(j);
    }
  }
  out.fill_diag(0);

  return out;
}

// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>

using namespace Rcpp;

//' @title Backward-in-time migration matrix conversion
//' @description Computes the backward-in-time migration matrix associated with a given forward-in-time migration matrix, coalescent rate and time scale
//' @param bit_mm NumericMatrix Forward-in-time migration matrix
//' @param coal_rate NumericVector Coalescent rates
//' @returns NumericMatrix Backward-in-time migration matrix
//' @export
//'
// [[Rcpp::export]]

NumericMatrix BitMigMatC(NumericMatrix fit_mm, NumericVector coal_rate) {
  int n_deme = fit_mm.nrow();
  NumericMatrix out(n_deme, n_deme);

  for (int i = 0; i < n_deme; ++i){
    for (int j = 0; j < n_deme; ++j){
      out(i,j) = fit_mm(j,i) * coal_rate(i) / coal_rate(j);
    }
  }
  out.fill_diag(0);

  return out;
}
