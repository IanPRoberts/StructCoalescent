// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace Rcpp;

//' Sample path from the distribution of an endpoint-conditioned CTMC
//'
//' Simulate a sample path from an endpoint conditioned CTMC by uniformization
//' using a pre-computed transition probability matrix.
//' Taken from ECctmc-0.2.5 under a GPL-3 license
//'
//' @param a,b States at the interval endpoints, provided as integers
//'    corresponding to rows of the CTMC rate matrix.
//' @param t0,t1 times of the interval endpoints
//' @param Q CTMC rate matrix
//' @param P CTMC transition probability matrix over the interval.
//'
//' @return matrix whose first column is the sequence of transition times
//' bookended by interval endpoints, and whose second column is the sequence of
//' states
//' @author Jon Fintzi (ECctmc-0.2.5)
// [[Rcpp::export]]

arma::mat sample_path(const int a, const int b, const double t0, const double t1, const arma::mat& Q, const arma::mat& P) {

  // Get the number of states and initialize vector of states
  int n_states = Q.n_rows;
  Rcpp::IntegerVector states = Rcpp::seq_len(n_states);

  // Get the length of the interval and the largest diagonal element of Q
  double T = t1 - t0;
  double m = max(abs(Q.diag()));

  // Construct the transition probability matrix and extract the a,b elem.
  double p_ab = P(a-1, b-1);

  // Generate the auxilliary transition matrix
  arma::mat R = arma::eye(n_states, n_states) + Q/m;

  // Sample threshold for determining the number of states
  Rcpp::NumericVector n_thresh = Rcpp::runif(1, 0, 1);

  // Initialize number of jumps and conditional probability of n jumps
  int n_jumps = 0;
  double c_prob = exp(-m*T) * (a == b) / p_ab;

  // proceed with sampling by uniformization
  // first the case when there are no jumps
  if(c_prob > n_thresh[0]) {

    // initialize matrix
    arma::mat path(2,2);

    // fill matrix
    path(0,0) = t0; path(0,1) = a;
    path(1,0) = t1;  path(1,1) = b;

    return path;

  } else {

    // increment the number of jumps and compute c_prob
    n_jumps += 1;
    c_prob += exp(-m*T) * pow(m*T, n_jumps) / Rcpp::internal::factorial(n_jumps) * R(a-1, b-1) / p_ab;

    // if there is exactly one jump
    if(c_prob > n_thresh[0]) {

      // if the endpoints match, the only jump is a virtual one
      if(a == b) {

        // initialize matrix
        arma::mat path(2,2);

        // fill matrix
        path(0,0) = t0; path(0,1) = a;
        path(1,0) = t1;  path(1,1) = b;

        return path;
        // if the endpoints don't match, the jump is real
                              } else {

                                        // initialize matrix
                                        arma::mat path(3,2);

                                        // fill matrix
                                        path(0,0) = t0; path(0,1) = a;
                                        path(1,0) = Rcpp::runif(1,t0,t1)[0]; path(1,1) = b;
                                        path(2,0) = t1; path(2,1) = b;

                                        return path;
                              }

                              // Else, there are at least two jumps
                    } else {

                              // Initialize a cube for storing powers of R
                              arma::cube R_pow(n_states, n_states, 8);
                              int R_pow_size = R_pow.n_slices;
                              R_pow.slice(0) = arma::eye(size(R));
                              R_pow.slice(1) = R;

                              // Initialize a vector for storing the transition probabilities
                              Rcpp::NumericVector state_probs(n_states);

                              // keep calculating the conditional probability of n jumps until
                              // the threshold is exceeded. store powers of R accordingly.
                              while(c_prob < n_thresh[0]) {

                                        // increment the number of jumps
                                        n_jumps += 1;

                                        // check whether to add additional slices to the R_pow cube
                                        if(n_jumps == R_pow_size) {
                                                  R_pow.insert_slices(R_pow.n_slices, 8);
                                                  R_pow_size = R_pow.n_slices;
                                        }

                                        // Add the new power of R to the cube and calculate c_prob
                                        R_pow.slice(n_jumps) = R_pow.slice(n_jumps - 1) * R;
                                        c_prob += exp(-m*T) * pow(m*T, n_jumps) / Rcpp::internal::factorial(n_jumps) * R_pow.slice(n_jumps)(a-1, b-1) / p_ab;
                              }

                              // initialize the path matrix
                              int path_nrows = n_jumps + 2;
                              arma::mat path(path_nrows, 2);
                              path(0,0) = t0;
                              path(0,1) = a;
                              path(path_nrows - 1, 0) = t1;
                              path(path_nrows - 1, 1) = b;

                              // transition times are uniformly distributed in the
                              // interval. Sample them, sort them, and place in path.
                              arma::colvec transitions = Rcpp::runif(n_jumps, t0, t1);
                              std::sort(transitions.begin(), transitions.end());
                              path(arma::span(1,n_jumps), 0) = transitions;

                              // Sample the states at the transition times
                              for(int j = 1; j < n_jumps + 1; ++j) {
                                        state_probs = arma::trans(R(path(j-1, 1) - 1, span::all)) % R_pow.slice(n_jumps-j)(span::all, b-1) / R_pow(path(j-1, 1)-1, b-1, n_jumps-j+1);
                                        path(j, 1) = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs)[0];
                              }

                              // Determine which transitions are virtual transitions
                              arma::vec keep_inds(path_nrows, arma::fill::ones);
                              for(int j = 1; j < n_jumps + 1; ++j) {
                                        if(path(j, 1) == path(j-1, 1)) {
                                                  keep_inds[j] = 0;
                                        }
                              }

                              // create a matrix for the complete path without virtual jumps
                              arma::mat path_comp = path.rows(arma::find(keep_inds == 1));

                              return path_comp;
                    }
          }
}
