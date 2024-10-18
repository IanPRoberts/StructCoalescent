#' Coalescent Node Pie Charts
#'
#' Tabulates the frequency of each deme at each leaf and coalescent node of a phylogenetic tree
#'
#'
#' @param strphylo_list List of dated structured phylogenetic trees on the same dated phylogenetic tree
#' @param plot_strphylo (optional) Dated phylogenetic tree on which to superimpose pie charts
#' @param plot (logical) Plot pie charts on top of plot_strphylo?
#' @param ... Additional parameters to pass to plot.strphylo
#'
#' @return Returns list giving strphylo corresponding to plotted structured phylogeny (if plot=TRUE) or the first element of strphylo_list (if plot=FALSE), and matrix of deme frequencies at coalescent events in increasing age order from the MRCA.
#'
#' @export

coalescent_node_pie_charts <- function(strphylo_list, plot_strphylo = NA, plot=TRUE, pie_cex=0.5, ...){
  ED_list <- lapply(strphylo_list, as.ED)
  if (any(is.na(plot_strphylo))){
    plot_strphylo <- strphylo_list[[1]]
  }
  plot_ED <- as.ED(plot_strphylo)
  n_trees <- length(ED_list)
  n_deme <- max(sapply(ED_list, function(x) max(x[,5])))

  #Strip migration events from tree
  topology <- ED_list[[1]]
  topology <- topology[(is.na(topology[,3])) | (!is.na(topology[,4])),]
  topology[,2:4] <- topology[,7:9]
  topology[!is.na(topology[,3]), 5] <- 0

  #Construct deme frequency at each coalescent node in ascending node age
  deme_freq <- matrix(0, nrow(topology), n_deme)

  for (tree_id in n_trees : 1){ #Loop in reverse tree order to leave node_order as order(topology[,6]) after final iteration
    ED <- ED_list[[tree_id]]
    ED <- ED[(is.na(ED[,3])) | (!is.na(ED[,4])),]
    ED[,2:4] <- ED[,7:9]

    node_order <- order(ED[,6]) #Store entries of deme_freq in ascending age (root = 0, newest leaf = max(ED[,6]))

    for (row_id in 1 : nrow(ED)){ #Loop robust only for non-simultaneous events, i.e. coalescent nodes
      deme_freq[row_id, ED[node_order[row_id], 5]] <- deme_freq[row_id, ED[node_order[row_id], 5]] + 1
    }
  }

  #Remove leaf deme frequencies from deme_freq (should be a.s. one deme throughout run)
  deme_freq <- deme_freq[!is.na(topology[node_order, 4]),]

  plot_ED_coal_nodes <- !is.na(plot_ED[,4])
  plot_ED_coal_node_order <- order(plot_ED[plot_ED_coal_nodes, 6])


  rownames(deme_freq) <- plot_ED[plot_ED_coal_nodes, 1][plot_ED_coal_node_order]

  if (plot){
    plot(plot_ED, ...)
    pie_plot <- rowSums(deme_freq == 0) < n_deme - 1 #Logical on whether 100% same deme observed (in which case no pie chart plotted!)
    ape::nodelabels(node = as.numeric(rownames(deme_freq))[pie_plot],
               pie = deme_freq[pie_plot,]/rowSums(deme_freq[pie_plot,]),
               cex = pie_cex)
  }

  return(list(strphylo = as.strphylo(plot_ED), node_freq = deme_freq))
}
