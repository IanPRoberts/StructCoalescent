#' Default StructCoalescent prior parameters
#'
#' Computes the prior parameters for the default prior specification for StructCoalescent
#'
#' @param strphylo strphylo object giving an initial structured phylogeny
#' @param n_deme Number of demes in the problem
#' @param M (optional) prior target number of migration events in a migration history. Default value places prior mass on maximum parsimony migration histories
#'
#' @return Updated extended data object with the proposal from the migration birth move
#'
#' @export

default_priors <- function(strphylo, n_deme, M=fitch(strphylo)$min_migs){
  ED <- as.ED(strphylo)
  n_leaf <- length(strphylo$tip.label)

  ED_DD <- DemeDecompC(ED, max(ED[,5]), NodeIndicesC(ED))
  k <- rowSums(ED_DD$k)
  coal_const <- sum(k * (k-1) * ED_DD$time.increments) / 2

  tree_length <- sum(strphylo$edge.length)

  prior_parameters <- c('cr_shape'=1,
                        'cr_rate'=coal_const / n_deme / (n_leaf - 1),
                        'mm_shape'=1,
                        'mm_rate'=(n_deme-1) * tree_length / M)
  return(prior_parameters)
}
