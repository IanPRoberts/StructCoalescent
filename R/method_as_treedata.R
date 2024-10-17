#' @export
as.treedata.ED <- function(ED){
  strphylo <- as.strphylo(ED)
  treedata <- as.treedata.strphylo(strphylo)
  return(treedata)
}

#' @export
as.treedata.strphylo <- function(strphylo){
  class(strphylo) <- 'phylo'
  treedata <- treeio::as.treedata(strphylo)
  treedata@data <- treeio::tibble(type = strphylo$node.deme, node = 1:length(strphylo$node.deme))
  return(treedata)
}
