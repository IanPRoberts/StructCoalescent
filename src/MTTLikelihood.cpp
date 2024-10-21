// ONLY WORKS FOR REAL EIGENVALUES AND EIGENVECTORS

// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List MTT_transition_kernelC(NumericMatrix ED, NumericMatrix bit_rates, NumericVector node_indices,
                            NumericVector eigen_vals, NumericMatrix eigen_vecs, NumericMatrix inverse_vecs){
  double log_like = 0;
  int n_deme=bit_rates.nrow();

  int parent_row;
  double time_inc;
  int current_deme;
  int parent_deme;
  double trans_prob;

  for (int i = 0; i < ED.nrow(); ++i){
    if (!NumericVector::is_na(ED(i,1))){ //If not root
      parent_row = node_indices[ED(i,1) - 1] - 1;
      time_inc = ED(i,5) - ED(parent_row, 5);

      current_deme = ED(i, 4) - 1;
      parent_deme = ED(parent_row, 4) - 1;

      log_like += bit_rates(current_deme, current_deme) * time_inc;

      if (parent_deme != current_deme){
        log_like += log(bit_rates(current_deme, parent_deme));
      }
    }

    if (NumericVector::is_na(ED(i, 2)) | !NumericVector::is_na(ED(i, 3))){ //If is_leaf | is_coalescent_node
      if (!NumericVector::is_na(ED(i,1))){ //If not root
        parent_row = node_indices[ED(i,6) - 1] - 1;
        time_inc = ED(i,5) - ED(parent_row, 5);

        current_deme = ED(i, 4) - 1;
        parent_deme = ED(parent_row, 4) - 1;

        trans_prob = 0;

        for (int a = 0; a < n_deme; ++a){
          trans_prob += exp(eigen_vals[a] * time_inc) * eigen_vecs(current_deme, a) * inverse_vecs(a, parent_deme);
        }

        log_like += -log(trans_prob);
      }
    }
  }

  List out = List::create(_["log.likelihood"] = log_like , _["likelihood"] = exp(log_like));
  return out;
}
