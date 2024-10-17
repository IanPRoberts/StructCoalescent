#' Evolutionary Parameter Gibbs Updates
#'
#' Gibbs updates for evolutionary parameters governing the structured coalescent
#'
#'
#' @param strphylo Dated structured phylogeny
#' @param prior_shape Gamma shape parameter for prior
#' @param prior_rate Gamma rate parameter for prior
#'
#' @return Updated parameters
#'
#' @export

coal_rate_gibbs_update <- function(strphylo, n_deme, prior_shape=1, prior_rate=1){
  return(ED_cr_gibbs(as.ED(strphylo), n_deme, NodeIndicesC(ED), prior_shape, prior_rate))
}

#' @rdname coal_rate_gibbs_update
#' @export

bit_mig_mat_gibbs_update <- function(strphylo, n_deme, prior_shape=1, prior_rate=1){
  return(ED_bmm_gibbs(as.ED(strphylo), n_deme, NodeIndicesC(ED), prior_shape, prior_rate))
}

ED_cr_gibbs <- function(ED, n_deme, ED_NI, prior_shape=1, prior_rate=1){
  c <- NodeCountC(ED, n_deme, ED_NI)$c
  DemeDecomp <- DemeDecompC(ED, n_deme, ED_NI)
  k <- DemeDecomp$k
  rate_consts <- colSums((k * (k-1) / 2) * DemeDecomp$time.increments)

  proposal <- rgamma(n_deme, prior_shape + c, prior_rate + rate_consts)
  return(proposal)
}


ED_bmm_gibbs <- function(ED, n_deme, ED_NI, prior_shape=1, prior_rate=1){
  m <- NodeCountC(ED, n_deme, ED_NI)$m
  deme_length <- sapply(1:n_deme, function(x){
    rows_in_deme <- ED[,5] == x
    sum(ED[rows_in_deme, 6] - ED[ED_NI[ED[rows_in_deme,2]], 6], na.rm=TRUE)
  })

  proposal <- matrix(rgamma(n_deme^2, prior_shape + m, prior_rate + deme_length), n_deme, n_deme)
  diag(proposal) <- 0
  return(proposal)
}
