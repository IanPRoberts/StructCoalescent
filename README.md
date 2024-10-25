
# StructCoalescent

<!-- badges: start -->
<!-- badges: end -->

The goal of StructCoalescent is to perform Bayesian inference using the
exact structured coalescent model.

## Installation

You can install StructCoalescent from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("IanPRoberts/StructCoalescent")
```

The package can be loaded using:

``` r
library(StructCoalescent)
```

## Example

This is a basic example of usage. First, we generate a heterochronous
dated structured phylogenetic tree with $n=10$ leaves sampled between
2015 and 2020 from $d=3$ demes.

``` r
n <- 10
d <- 3
leaf_data <- cbind('Leaf_ID'=paste0('L', 1:n),
                   'Leaf_age'=sample(2018:2020, n, TRUE),
                   'Leaf_deme'=sample(1:d, n, TRUE))
coalescent_rates <- rexp(d, 1)
migration_matrix <- matrix(rexp(d^2, 20), d, d); diag(migration_matrix) <- 0


strphylo <- rstrphylo(n, d, coalescent_rates, migration_matrix, leaf_data)
plot(strphylo, time_axis=TRUE, root_time=strphylo$root.age, show.tip.label=TRUE)
```

<img src="man/figures/README-rstrphylo-1.png" width="100%" />

We can run $N=100,000$ MCMC algorithm taking this dated structured
phylogeny as input with default prior distributions as follows (NOT
RUN):

``` r
StructCoalescent_mcmc(N=1e4, strphylo, coalescent_rates, migration_matrix, stdout_log=FALSE, thin=100, proposal_rates=c(20, 1, 1))
```

MCMC output consists of three files (by default added to current working
directory `getwd()`):

- `StructCoalescent.trees` contains a thinned MCMC sample of dated
  structured phylogenies in a BEAST meta-commented Newick format, which
  can be read with `treeio::read.beast`.
- `StructCoalescent.log` contains a thinned MCMC sample of evolutionary
  parameters, likelihood evaluations and the current radius of the
  subtree in a csv format, which can be read with `read.csv`
- `StructCoalescent.freq` contains details of the total number of
  accepted proposals of each type as well as the total number of
  attempted proposals of that type.

We can present a sample of migration histories using a consensus
migration history, with sections of the branch assigned a deme provided
a proportion $p$ of sampled migration histories observe the same deme at
that position as follows:

``` r
strphylo_list <- lapply(treeio::read.beast('./StructCoalescent.trees'), as.strphylo)
consensus_strphylo <- exact_consensus(strphylo_list, consensus_prob=0.75)
```

<img src="man/figures/README-consensus-history-1.png" width="100%" />

Evolutionary parameter samples can be presented in a variety of ways
including trace plots and histograms:

### Trace plots

``` r
StructCoalescent_log <- read.csv('./StructCoalescent.log',
                                 header=TRUE)
layout(matrix(1:d^2, d, d))
for (col_id in 1:d){
  for (row_id in 1:d){
    if (col_id == row_id){
      col_name <- paste0('coal_rate_', col_id)
      title <- list(paste0('theta[', col_id, ']'))
      true_param <- coalescent_rates[col_id]
    } else{
      col_name <- paste('backward_migration_rate', row_id, col_id, sep='_')
      title <- paste0('lambda[', row_id, ',', col_id, ']')
      true_param <- migration_matrix[row_id, col_id]
    }
    plot(StructCoalescent_log$sample, StructCoalescent_log[,col_name], type='l', xlab='Sample', ylab=title)
    abline(h=true_param, lty=2, col='red')
  }
}
```

<img src="man/figures/README-evolutionary-parameter-traces-1.png" width="100%" />

### Histograms

``` r
StructCoalescent_log <- read.csv('./StructCoalescent.log',
                                 header=TRUE)
layout(matrix(1:d^2, d, d))
for (col_id in 1:d){
  for (row_id in 1:d){
    if (col_id == row_id){
      col_name <- paste0('coal_rate_', col_id)
      title <- list(paste0('theta[', col_id, ']'))
      true_param <- coalescent_rates[col_id]
    } else{
      col_name <- paste('backward_migration_rate', row_id, col_id, sep='_')
      title <- paste0('lambda[', row_id, ',', col_id, ']')
      true_param <- migration_matrix[row_id, col_id]
    }
    hist(StructCoalescent_log[,col_name], xlab=title, freq=FALSE, main='')
    abline(v=true_param, lty=2, col='red')
  }
}
```

<img src="man/figures/README-evolutionary-parameter-hists-1.png" width="100%" />
