// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
using namespace Rcpp;

//' @title NodeCountC
//' @description Counts the number of migrations between pairs of demes, and coalescences within each deme
//' @param ED NumericMatrix Extended data object representing structured phylogeny
//' @param n_deme int Number of demes modelled under ED
//' @param node_indices NumericVector Vector giving row numbers in ED for node labels
//' @returns List containing matrix m with element (i,j) giving the number of migrations i -> j backwards in time, and a vector c with element i giving the number of coalescences occurring in deme i
//'
//' @export
//'
// [[Rcpp::export]]

List NodeCountC(NumericMatrix ED, int n_deme, NumericVector node_indices) {
  NumericMatrix m(n_deme, n_deme);
  NumericVector c(n_deme);

  int nrow = ED.nrow();
  int child_row = 0;
  int mig_target = 0;
  int mig_origin = 0;

  //Could start at i = root_row to skip all leaf checks

  for (int i = 0; i < nrow; ++i) { //Find root row
    if(NumericVector::is_na(ED(i,2))){ //Leaf
      continue;
    } else if (NumericMatrix::is_na(ED(i, 3)) == 0){ //Coalescence
      c[ED(node_indices[i] - 1, 4) - 1] += 1;
    } else { //Migration
      child_row = node_indices[ED(i, 2) - 1] - 1;
      mig_target = ED(i, 4) - 1;
      mig_origin = ED(child_row, 4) - 1;
      m(mig_origin, mig_target) += 1;
    }
  }

  List out = List::create(_["m"] = m, _["c"] = c);

  return out;
}
