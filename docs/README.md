---
title: "README"
output: html_document
---

# StructCoalescent

<!-- badges: start -->
<!-- badges: end -->

StructCoalescent is an R package designed to perform Bayesian inference using the exact structured coalescent model.
A paper describing the methods behind StructCoalescent has been published, see [Roberts I, Everitt RG, Koskela J, Didelot X (2025) Bayesian Inference of Pathogen Phylogeography using the Structured Coalescent Model. PLOS Computational Biology 21(4): e1012995](https://doi.org/10.1371/journal.pcbi.1012995).


## Installation

StructCoalescent can be installed from [GitHub](https://github.com/IanPRoberts/StructCoalescent) using devtools using:

```{r gh-installation, eval = FALSE}
if ( !require( devtools, quietly = TRUE ) ){
  install.packages("devtools")
}
devtools::install_github("IanPRoberts/StructCoalescent")
```

and can then be loaded using

```{r}
library(StructCoalescent)
```

## Examples

Further examples and documentation on how to use StructCoalescent:

1. [Converting phylogenetic trees](Converting_file_types.html)
2. [Setting up MCMC](Setting_up_MCMC.html)
3. Summarising MCMC results
