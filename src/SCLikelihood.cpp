// [[Rcpp::interfaces(r, cpp)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

double SC_like_C(NumericMatrix ED, NumericVector coal_rate, NumericMatrix bit_mig_mat, NumericVector node_indices) {
  double n_deme = coal_rate.size();
  double like = 0;
  int row;
  int current_deme;
  int child_deme;
  double time_inc;
  NumericVector k(n_deme); //Contemporary lineages
  NumericVector mm_rowsum = rowSums(bit_mig_mat);

  arma::vec event_times = ED.column(5);
  arma::uvec node_order = arma::sort_index(event_times, "descend"); //order(ED[,6], decreasing = TRUE) in R

  for (int id = 0; id < ED.nrow() - 1; ++id){
    row = node_order[id];
    current_deme = ED(row, 4) - 1;

    if (NumericVector::is_na(ED(row, 2))){ //Leaf
      k(current_deme) += 1;
    } else if (NumericVector::is_na(ED(row, 3))){ //Migration
      child_deme = ED(node_indices[ED(row, 2) - 1] - 1, 4) - 1;
      k(child_deme) -= 1;
      k(current_deme) += 1;

      like += log(bit_mig_mat(child_deme, current_deme)); //Normalising constant
    } else { // Coalescence
      k(current_deme) -= 1;
      like += log(coal_rate[current_deme]); //Normalising constant
    }

    time_inc = ED(row, 5) - ED(node_order[id + 1], 5); //Increment until following event (mig/coal/leaf)
    like -= time_inc * (sum(k * (k - 1) * coal_rate) / 2 + sum(k * mm_rowsum));
  }

  like += log(coal_rate[current_deme]); //Normalising constant for root coalescence (const added when node is child of edge)

  return(like);
}
