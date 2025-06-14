#' @export

plot.strphylo <- function (strphylo, n_deme = NA, time_axis = FALSE, root_time = NA, color.palette=NULL,
                           ...){
  edge <- strphylo$edge
  if (is.na(n_deme)) {
    n_deme <- max(strphylo$node.deme)
  }

  # if (n_deme <= 7){
  #   color.palette <- c(palette.colors(n_deme + 1, palette = 'R4')[-1], "black")
  # } else {
  #   color.palette <- c(palette.colors(n_deme, palette = 'Polychrome 36'), "black")
  # }
  if (is.null(color.palette)){
    color.palette <- c(rainbow(n_deme), 'black')
  }
  strphylo$node.deme[strphylo$node.deme == 0] <- n_deme + 1
  edge.color <- color.palette[strphylo$node.deme[edge[, 2]]]
  ape::plot.phylo(strphylo, edge.color = edge.color,
       ...)
  if (time_axis) {
    if (is.na(root_time)) {
      root_time <- 0
    }
    ape::axisPhylo(root.time = root_time, backward = FALSE)
  }
}

#' @export
plot.ED <- function(ED, n_deme=NA, time_axis=FALSE, root_time=NA, ...){
  strphylo <- as.strphylo(ED)
  plot.strphylo(strphylo, n_deme, time_axis, root_time, ...)
}
