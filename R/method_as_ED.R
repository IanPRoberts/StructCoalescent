as.ED <- function (x, ...)
{
  if (identical(class(x), "ED")) return(x)
  UseMethod("as.ED")
}

#' @export
as.ED.default <- function(x, ...){
  if (inherits(x, "ED")) return(x)
  stop('Object does not inherit the class "ED": found no appropriate method to convert it')
}

#' @export
as.ED.strphylo <- function(strphylo){
  n <- 1 + nrow(strphylo$edge) #Number of lineages
  ED <- cbind('Node_ID'=1:n,
              'Parent'=NA,
              'Child_1'=NA,
              'Child_2'=NA,
              'Deme'=strphylo$node.deme,
              'Node_Age'=ape::node.depth.edgelength(strphylo),
              'Parent_coal'=NA,
              'Child_coal_1'=NA,
              'Child_coal_2'=NA)

  for (i in 1 : nrow(strphylo$edge)){
    if (ED[strphylo$edge[i,1],6] < ED[strphylo$edge[i,2],6]){
      min_age <- 1
      max_age <- 2
    } else {
      min_age <- 2
      max_age <- 1
    }
    ED[strphylo$edge[i, max_age], 2] <- strphylo$edge[i, min_age]

    if (is.na(ED[strphylo$edge[i, min_age],3])){
      ED[strphylo$edge[i, min_age],3] <- strphylo$edge[i, max_age]
    } else {
      ED[strphylo$edge[i, min_age],4] <- strphylo$edge[i, max_age]
    }
  }

  ED[,7:9] <- ED[,2:4]
  #Parent_coal
  while (sum(is.na(ED[ED[,7], 4])) > 1){
    ED[,7][is.na(ED[ED[,7], 4])] <- ED[ED[is.na(ED[ED[,7], 4]), 7], 2]
  }

  #Child_coal_1
  active_rows <- is.na(ED[ED[,8],4]) & !is.na(ED[ED[,8],3])
  while(sum(active_rows) > 0){
    ED[active_rows, 8] <- ED[ED[active_rows, 8], 3]
    active_rows <- is.na(ED[ED[,8],4]) & !is.na(ED[ED[,8],3])
  }

  #Child_coal_2
  active_rows <- is.na(ED[ED[,9],4]) & !is.na(ED[ED[,9],3])
  while(sum(active_rows) > 0){
    ED[active_rows, 9] <- ED[ED[active_rows, 9], 3]
    active_rows <- is.na(ED[ED[,9],4]) & !is.na(ED[ED[,9],3])
  }

  # Redo labels for leaf 1:n, root n+1, coal nodes (n+2) : (2n-1), mig nodes > 2n-1
  leaf_rows <- is.na(ED[,3])
  coalescence_rows <- !is.na(ED[,4])
  migration_rows <- !(leaf_rows | coalescence_rows)

  ED <- rbind(ED[leaf_rows,],
              ED[coalescence_rows,],
              ED[migration_rows,])

  max_label <- max(ED[,1], na.rm=TRUE)
  ED[,1:4] <- ED[,1:4] + max_label
  for (i in 1 : n){#nrow(ED)){
    ED[,1:4][ED[,1:4] == ED[i,1]] <- i
  }
  class(ED) <- "ED"
  return(ED)
}


#' @export
as.ED.phylo <- function(phylo, node.deme){
  strphylo <- as.strphylo(phylo, node.deme)
  ED <- as.ED.strphylo(strphylo)
  return(ED)
}

#' @export
as.ED.matrix <- function(x){
  if (ncol(x) == 9){
    class(x) <- 'ED'
    return(x)
  }
  stop('Object does not inherit the class "ED": found no appropriate method to convert it')
}
