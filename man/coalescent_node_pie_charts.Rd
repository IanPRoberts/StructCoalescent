% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coal_node_pie_charts.R
\name{coalescent_node_pie_charts}
\alias{coalescent_node_pie_charts}
\title{Coalescent Node Pie Charts}
\usage{
coalescent_node_pie_charts(
  strphylo_list,
  plot_strphylo = NA,
  plot = TRUE,
  pie_cex = 0.5,
  ...
)
}
\arguments{
\item{strphylo_list}{List of dated structured phylogenetic trees on the same dated phylogenetic tree}

\item{plot_strphylo}{(optional) Dated phylogenetic tree on which to superimpose pie charts}

\item{plot}{(logical) Plot pie charts on top of plot_strphylo?}

\item{...}{Additional parameters to pass to plot.strphylo}
}
\value{
Returns list giving strphylo corresponding to plotted structured phylogeny (if plot=TRUE) or the first element of strphylo_list (if plot=FALSE), and matrix of deme frequencies at coalescent events in increasing age order from the MRCA.
}
\description{
Tabulates the frequency of each deme at each leaf and coalescent node of a phylogenetic tree
}
