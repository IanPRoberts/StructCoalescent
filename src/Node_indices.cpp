// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>

using namespace Rcpp;

//' @title Node Indices
//' @description Returns a vector giving the row indices of the labels in ED
//' @param ED NumericMatrix Extended data object representing structured coalescent genealogy
//'
//' @export
//'
// [[Rcpp::export]]

NumericVector NodeIndicesC(NumericMatrix ED) {
  int max_lab = max(ED(_,0));

  NumericVector node_indices(max_lab);

  for (int i = 0; i < ED.rows(); ++i){
    node_indices[ED(i,0) - 1] = i+1;
  }
  return node_indices;
}
