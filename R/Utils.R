#' The Truncated Exponential Distribution
#'
#' Density, distribution function and random generation for the exponential
#' distribution with rate \code{rate} truncated above by \code{c}
#'
#' @param x,q vector of quantiles
#' @param n number of observations
#' @param rate vector of rates
#' @param c truncation cutoff
#'
#' @export

#Density function
dexp_trunc <- function(x, c = Inf, rate = 1) rate * exp(-rate * x) / (1 - exp(-rate * c))

#' @rdname dexp_trunc
#' @export
#Distribution function
pexp_trunc <- function(q, c = Inf, rate = 1) (1 - exp(-rate * q)) / (1 - exp(-rate*c))

#' @rdname dexp_trunc
#' @export
#i.i.d. sample of size n
rexp_trunc <- function(n, c = Inf, rate = 1) - (1/rate) * log(1 - (1 - exp(-rate * c)) * runif(n))
