#' @export

plot.strphylo <- function (strphylo, n_deme = NA, time_axis = FALSE, root_time = NA,
                           ...){
  edge <- strphylo$edge
  if (is.na(n_deme)) {
    n_deme <- max(strphylo$node.deme)
  }
  # color.palette <- c(palette.colors(n_deme, palette = 'Polychrome 36'), "black")
  color.palette <- c(rainbow(n_deme), 'black')
  strphylo$node.deme[strphylo$node.deme == 0] <- n_deme + 1
  edge.color <- color.palette[strphylo$node.deme[edge[, 2]]]
  # class(strphylo) <- 'phylo'
  plot.phylo(strphylo, edge.color = edge.color, edge.width = 2, show.tip.label = FALSE,
       ...)
  if (time_axis) {
    if (is.na(root_time)) {
      root_time <- 0
    }
    axisPhylo(root.time = root_time, backward = FALSE)
  }
}
