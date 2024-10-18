#' @title Structured Coalescent likelihood
#' @description Computes the likelihood of a structured phylogenetic tree under the structured coalescent model
#' @param strphylo Structured dated phylogeny
#' @param bit_mig_mat Matrix of backwards-in-time migration rates between pairs of demes
#' @param coal_rate Vector of coalescent rate in each deme (1 / effective population size)
#'
#' @returns List containing log-likelihood and likelihood of the structured phylogenetic tree
#' @export

structured_coalescent_likelihood <- function(strphylo, coal_rate, bit_mig_mat){
  ED <- as.ED(strphylo)
  diag(bit_mig_mat) <- 0
  log_likelihood <- SC_like_C(ED, coal_rate, bit_mig_mat, NodeIndicesC(ED))

  return(list(likelihood=exp(log_likelihood), log_likelihood=log_likelihood))
}

#' @title DTA likelihood
#' @description Computes the likelihood of a structured phylogenetic tree under the Discrete Trait Analysis (DTA) model
#' @param strphylo Structured dated phylogeny
#' @param fit_mig_mat Matrix of forwards-in-time migration rates between pairs of demes
#' @param bit_mig_mat Matrix of backwards-in-time migration rates between pairs of demes
#' @param coal_rate Vector of coalescent rate in each deme (1 / effective population size)
#'
#' @returns List containing log-likelihood and likelihood of the structured phylogenetic tree
#' @export

DTA_likelihood <- function(strphylo, fit_mig_mat=FitMigMatC(bit_mig_mat, coal_rate), bit_mig_mat=BitMigMatC(fit_mig_mat, coal_rate), coal_rate){
  ED <- as.ED(strphylo)
  diag(fit_mig_mat) <- 0
  likelihood <- DTALikelihoodC(ED, fit_mig_mat, NodeIndicesC(ED))

  return(list(likelihood=likelihood$likelihood, log_likelihood=likelihood$log.likelihood))
}
