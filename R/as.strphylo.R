as.strphylo <- function (x, ...){
  if (identical(class(x), "strphylo")) return(x)
  UseMethod("as.strphylo")
}

#' @export
as.strphylo.default <- function(x, ...){
  if (inherits(x, "strphylo")) return(x)
  stop('Object does not inherit the class "strphylo": found no appropriate method to convert it')
}

#' @export
as.strphylo.ED <- function(ED){
  n.nodes <- dim(ED)[1]
  n.tips <- sum(is.na(ED[,3]))

  #Check no gaps in node labelling scheme to allow phylo object to be generated
  if (max(ED[,1]) > n.nodes){
    missing.labels <- (1:n.nodes)[! (1:n.nodes) %in% ED[,1]]
    extra.labels <- as.vector(ED[ED[,1] > n.nodes,1])  #as.vector needed in case only 1 extra label has appeared

    node.label.mat <- ED[,1:4]
    count <- 1
    for (i in extra.labels){
      node.label.mat[node.label.mat == i] <- missing.labels[count]
      count <- count + 1
    }
    ED[, 1:4] <- node.label.mat
  }
  edge.list <- list()
  edge.length <- numeric(0)
  count <- 1
  for (i in (1:n.nodes)[-(n.tips + 1)]){
    row <- which(ED[,1] == i)
    parent.row <- which(ED[,1] == ED[row, 2])
    edge.list[[count]] <- c(ED[row, 2], i)
    edge.length[count] <- ED[row, 6] - ED[parent.row, 6]
    count <- count + 1
  }
  edge <- do.call(rbind,edge.list)  #construct edge matrix from edge.list

  phylo <- list()
  class(phylo) <- c('strphylo', 'phylo')
  phylo$edge <- edge
  phylo$edge.length <- edge.length
  phylo$tip.label <- 1:n.tips
  phylo$Nnode <- n.nodes - n.tips
  phylo$node.deme <- ED[order(ED[,1]),5]  #Order supplies the ordering of the rows in ED to get node demes in correct order

  phylo <- ape::ladderize(phylo)

  return(phylo)
}

#' @export
as.strphylo.phylo <- function(phylo, node.deme){
  phylo$node.deme <- node.deme
  class(phylo) <- c('strphylo', 'phylo')
  return(phylo)
}
