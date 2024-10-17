// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>

using namespace Rcpp;

//' @title DTALikelihoodC
//' @description Computes the DTA likelihood of a structured genealogy
//' @param ED NumericMatrix Extended data object representing structured phylogeny
//' @param coal_rate NumericVector Coalescent rate in each deme (1 / effective population size)
//' @param bit_mig_mat NumericMatrix Backwards-in-time migration rates between pairs of demes
//' @param node_indices NumericVector Vector giving row numbers in ED for node labels
//' @returns List containing log-likelihood and likelihood of the structured coalescent genealogy
//'
//' @export
//'
// [[Rcpp::export]]

List DTALikelihoodC(NumericMatrix ED, NumericMatrix fit_mig_mat, NumericVector node_indices){
  double log_like = 0;
  int n_deme = fit_mig_mat.nrow();

  NumericVector fit_mm_rowsum(n_deme);
  for (int i = 0; i < n_deme; ++i){
    fit_mm_rowsum(i) = sum(fit_mig_mat(i, _));
  }

  int parent_row;
  double time_inc;
  int parent_deme;
  int node_deme;

  for (int i = 0; i < ED.nrow(); ++i){
    if (!NumericVector::is_na(ED(i,1))){
      parent_row = node_indices[ED(i,1) - 1] - 1;
      time_inc = ED(i,5) - ED(parent_row, 5);

      parent_deme = ED(i, 4) - 1;

      if (NumericVector::is_na(ED(i,2))){ //Leaf node
        node_deme = parent_deme;
      } else {
        node_deme = ED(node_indices(ED(i,2) - 1) - 1, 4) - 1;
      }

      log_like -= fit_mm_rowsum[parent_deme] * time_inc;

      if (parent_deme != node_deme){
        log_like += log(fit_mig_mat(parent_deme, node_deme));
      }
    }
  }

  List out = List::create(_["log.likelihood"] = log_like , _["likelihood"] = exp(log_like));
  return out;
}
