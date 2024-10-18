#' Fitch Algorithm
#'
#' Computes the number of required migration events for a maximum parsimony migration history on a given structured phylogenetic tree
#'
#' @param strphylo Dated structured phylogenetic tree containing at least tip demes
#'
#' @export

fitch <- function(strphylo){
  ED <- as.ED(strphylo)
  fitch_min <- fitch_ED(ED)
  return(list(min_migs = fitch_min$min_migs, strphylo = as.strphylo(fitch_min$ED)))
}

fitch_ED <- function(ED){
  n_deme <- max(ED[,5]) #Maximum number of required demes = max observed deme

  topology <- ED[is.na(ED[,3]) | !is.na(ED[,4]), ]
  topology[, 2:4] <- topology[, 7:9]
  top_NI <- NodeIndicesC(topology)

  optimal_demes <- matrix(FALSE, nrow(topology), n_deme)
  active_rows <- top_NI[topology[is.na(topology[,3]), 1]] #Start with leaf nodes as active rows

  ####### Root-ward sweep
  while(length(active_rows) > 0){
    for (row in active_rows){
      if (is.na(topology[row, 3])){ #Leaf optimal state == sampled deme
        optimal_demes[row, topology[row, 5]] <- TRUE
      } else { #Coalescent node optimal states either intersection (if not empty) or union of child optimal states
        child_rows <- top_NI[topology[row, 3:4]]
        intersection <- optimal_demes[child_rows[1],] & optimal_demes[child_rows[2],]
        if (any(intersection)){ #Intersection not empty
          optimal_demes[row,] <- intersection
        } else {
          optimal_demes[row,] <- optimal_demes[child_rows[1],] | optimal_demes[child_rows[2],]
        }
      }
    }
    active_rows <- na.omit(top_NI[topology[active_rows, 2]])
  }

  ######## Leaf-ward sweep
  root_row <- top_NI[topology[is.na(topology[,2]), 1]]
  topology[root_row, 5] <- sample.int(n_deme, 1, prob=optimal_demes[root_row,])

  min_migs <- 0
  active_rows <- topology[root_row, 3:4]

  while (length(active_rows) > 0){
    for (row in active_rows){
      parent_row <- top_NI[topology[row, 2]]
      parent_deme <- topology[parent_row, 5]

      if (optimal_demes[row, parent_deme]){
        topology[row, 5] <- parent_deme
      } else {
        topology[row, 5] <- sample.int(n_deme, 1, prob=optimal_demes[row,])
        min_migs <- min_migs + 1
      }
    }
    active_rows <- na.omit(top_NI[topology[active_rows, 3:4]])
  }

  return(list(min_migs = min_migs, ED = structure(topology, class='ED')))
}
