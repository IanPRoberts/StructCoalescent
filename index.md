# StructCoalescent

<!-- badges: start -->
<!-- badges: end -->
  
StructCoalescent is an R package designed to perform Bayesian inference using the exact structured coalescent model.
A paper describing the methods behind StructCoalescent has been published, see [Roberts I, Everitt RG, Koskela J, Didelot X (2025) Bayesian Inference of Pathogen Phylogeography using the Structured Coalescent Model. PLOS Computational Biology 21(4): e1012995](https://doi.org/10.1371/journal.pcbi.1012995).


## Installation

StructCoalescent can be installed from [GitHub](https://github.com/IanPRoberts/StructCoalescent) using:
  
  ``` r
if ( !require( pak, quietly = TRUE ) ){
  install.packages("pak")
}
pak::pak("IanPRoberts/StructCoalescent")
```

and can then be loaded using

``` r
library(StructCoalescent)
```