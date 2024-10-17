#' Majority Consensus Migration History
#'
#' Computes a majority consensus migration history for a sample of structured phylogenies on the same underlying phylogeny
#'
#'
#' @param strphylo_list Sample of migration histories for the same dated phylogeny
#' @param consensus_prob Consensus probability (must be greater than 0.5 to guarantee unique consensus history)
#' @param plot Logical value whether to output plot of consensus migration history
#' @param ... Additional arguments to pass to plot
#'
#' @return Consensus migration history in ED format
#'
#' @export
exact_consensus <- function(strphylo_list, consensus_prob=0.5, plot=TRUE, ...){
  ED_list <- lapply(strphylo_list, as.ED)
  con_ED <- exact_consensus_ED(ED_list, consensus_prob)

  if (plot){
    plot(con_ED, ...)
  }
  return(as.strphylo(con_ED))
}


exact_consensus_ED <- function(ED_list, consensus_prob=0.5){
  n_deme <- max(sapply(ED_list, function(x) max(x[,5])))
  n_trees <- length(ED_list)
  NI_list <- lapply(ED_list, NodeIndicesC)

  if (consensus_prob < 0.5){
    stop("Consensus probability won't give unique consensus tree")
  } else if (consensus_prob == 0.5){
    consensus_prob <- 0.501
  }

  consensus_freq <- n_trees * consensus_prob

  coal_times <- sort(unique(ED_list[[1]][!is.na(ED_list[[1]][,4]), 6]))

  con_ED <- ED_list[[1]]
  con_ED[, 2:4] <- con_ED[, 7:9]
  con_ED <- con_ED[(!is.na(con_ED[,4])) | is.na(con_ED[,3]), ]
  con_ED[!is.na(con_ED[,3]),5] <- 0
  con_NI <- NodeIndicesC(con_ED)

  new_label <- max(con_ED[,1]) + 1

  leaf_labels <- con_ED[is.na(con_ED[,3]), 1]
  n_leaf <- length(leaf_labels)

  coal_events <- !is.na(con_ED[,4])
  coal_labels <- con_ED[coal_events,c(1,6)][order(con_ED[coal_events, 6], decreasing = TRUE),] #Sort coal labels into descending


  current_coal_labels <- NA

  for (label_id in 1 : (2 * n_leaf - 2)){ # Iterate to 2n-2 to not try and look at edge above root
    if (label_id <= n_leaf){
      current_rows <- sapply(NI_list, function(x) x[leaf_labels[label_id]])
      child_coal <- leaf_labels[label_id]
    } else {
      coal_id <- label_id - n_leaf
      #Identify coalescent node labels using a.s. unique ordering in case of label switching and numerical precision
      current_coal_labels <- sapply(ED_list, function(ED){
        coal_events <- !is.na(ED[,4])
        coal_order <- order(ED[coal_events, 6], decreasing=TRUE)
        return(ED[coal_events, 1][coal_order][coal_id])
      })
      current_rows <- sapply(1:n_trees, function(x) NI_list[[x]][current_coal_labels[x]])
      child_coal <- coal_labels[label_id - n_leaf, 1]
    }

    con_ED_row <- con_NI[child_coal]
    current_demes <- sapply(1:n_trees, function(x){ED_list[[x]][current_rows[x], 5]})

    deme_freq <- sapply(1:n_deme, function(x) sum(current_demes == x))
    if (any(deme_freq >= consensus_freq)){
      con_ED[con_ED_row, 5] <- which(deme_freq >= consensus_freq)
      current_consensus_deme <- con_ED[con_ED_row, 5] # Current deme in consensus tree
    } else {
      current_consensus_deme <- 0 # Current deme in consensus tree
    }

    event_times <- do.call(rbind, lapply(1:n_trees, function(x){
      ED <- ED_list[[x]]
      #Migrations on branch have same parent coal as current row and child coal 1 equal to current row
      branch_migs <- (ED[,7] == ED[current_rows[x], 7]) & (ED[,8] == ED[current_rows[x], 1])
      branch_migs[is.na(branch_migs)] <- FALSE

      event_times <- matrix(FALSE, sum(branch_migs), n_trees + 1)
      event_times[, x+1] <- TRUE
      event_times[, 1] <- ED[branch_migs, 6]

      return(event_times)
    }))


    if (nrow(event_times) > 0){
      unique_times <- unique(event_times[,1])
      n_unique <- length(unique_times)
      if (n_unique < nrow(event_times)){
        #Combine duplicated events between trees into a single row
        new_event_times <- matrix(NA, n_unique, n_trees + 1)
        for (unique_id in 1 : length(unique_times)){
          non_unique_rows <- event_times[,1] == unique_times[unique_id]
          if (sum(non_unique_rows) == 1){
            new_event_times[unique_id,] <- event_times[non_unique_rows, ]
          } else {
            new_event_times[unique_id,] <- c(unique_times[unique_id],
                                             colSums(event_times[non_unique_rows, -1]))
          }
        }
        event_times <- new_event_times
      }

      if (nrow(event_times) > 1){
        event_times <- event_times[order(event_times[,1], decreasing = TRUE),] #Sort into decreasing event age
      }

      con_ED[con_ED_row, 2] <- new_label #Update parent of leaf node
      new_migs <- 0

      for (event_id in 1 : nrow(event_times)){
        current_ED <- as.logical(event_times[event_id, -1]) #Current trees to iterate event row over
        current_rows[current_ED] <- sapply((1:n_trees)[current_ED], function(x) NI_list[[x]][ED_list[[x]][current_rows[x], 2]])
        current_demes[current_ED] <- sapply((1:n_trees)[current_ED], function(x){ED_list[[x]][current_rows[x], 5]})

        deme_freq <- sapply(1:n_deme, function(x) sum(current_demes == x))
        if (any(deme_freq >= consensus_freq)){
          consensus_deme <- which(deme_freq >= consensus_freq)
        } else {
          consensus_deme <- 0
        }

        if (consensus_deme != current_consensus_deme){
          # Add new migration event if consensus deme different to current deme in consensus tree
          new_migs <- new_migs + 1
          con_ED <- rbind(con_ED,
                          c(new_label, #Node Id
                            new_label + 1, #Parent
                            new_label - 1, #con_ED[con_ED_row, 1], #Child 1
                            NA, # Child 2
                            consensus_deme, #Deme
                            event_times[event_id, 1], #Node Age
                            con_ED[con_ED_row, 7], #Parent coal
                            child_coal, #child coal 1
                            NA #child coal 2
                          ))
          con_NI[new_label] <- nrow(con_ED)
          new_label <- new_label + 1
        }

        current_consensus_deme <- consensus_deme
      }

      parent_coal_row <- con_NI[con_ED[con_NI[child_coal], 7]]
      which_child <- which(con_ED[parent_coal_row, 8:9] == child_coal)

      if (new_migs > 0){
        con_ED[nrow(con_ED), 2] <- con_ED[nrow(con_ED), 7] #Parent of final migration is coalescent node above
        con_ED[parent_coal_row, 2+which_child] <- new_label - 1 #Child of parent coalescent node

        con_ED[nrow(con_ED) + 1 - new_migs, 3] <- child_coal
      } else {
        con_ED[con_NI[child_coal], 2] <- con_ED[parent_coal_row, 1]
        con_ED[parent_coal_row, 2+which_child] <- child_coal
      }
    }
  }
  return(structure(con_ED, class='ED'))
}
