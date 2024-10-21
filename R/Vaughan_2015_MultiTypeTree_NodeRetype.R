#' MultiTypeTree Subtree
#'
#' Extracts a subtree (and associated sub-migration history) of a structured
#' genealogy consisting of the parent branch of a coalescent event and a fixed
#' number of generations of descendant branches
#'
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param st_depth Number of generations below the selected coalescent node to include in the subtree
#' @param node_indices Vector of row indices corresponding to which row of ED corresponds to each node label
#' @param selected_node (optional) Selected coalescent node to use as centre of subtree
#'
#' @return List consisting of ED, the structured phylogeny input, and st_labels, a reduced extended data structure consisting of the coalescent nodes in the subtree only
#'
#' @export

MTT_st_coal_node <- function(ED, st_depth = 1, node_indices = NodeIndicesC(ED), selected_node = NA){
  if (is.na(selected_node)){
    coal_nodes <- ED[!is.na(ED[,4]), 1]
    selected_node <- sample(coal_nodes, 1)
  }

  selected_row <- node_indices[selected_node]
  st_root <- ED[selected_row, 7]
  max_label <- max(ED[,1]) + 1

  if (is.na(st_root)){
    st_root <- selected_node
    st_labels <- st_root
  } else {
    st_labels <- c(st_root, selected_node)
  }

  active_rows <- selected_row

  for (gen in 1 : st_depth){
    child_coals <- unname(na.omit(ED[active_rows, 8:9]))
    st_labels <- append(st_labels,
                        child_coals)
    active_rows <- node_indices[child_coals]
  }


  return(list(ED = ED, st_labels = ED[node_indices[st_labels],]))
}

#' MultiTypeTree Node Retype
#'
#' Generates a new migration history on a subtree drawn using MTT_st_coal_node()
#'
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to coalescent nodes in the selected subtree only
#' @param bit_rates Transition matrix of a Markov process corresponding to a process with rates given by the backwards-in-time migration matrix
#' @param node_indices Vector of row indices corresponding to which row of ED corresponds to each node label
#'
#' @return  List consisting of proposal, a structured phylogeny with updates made to the selected subtree, and prop_prob, the probability of the given update

MTT_node_retype <- function(ED, st_labels, bit_rates, node_indices = NodeIndicesC(ED), eigen_vals = eigen(bit_rates)$values, eigen_vecs = eigen(bit_rates)$vectors, inverse_vecs = solve(eigen_vecs)){
  st_root <- unname(st_labels[1,1]) #st_labels[1,] always corresponds to st_root by construction
  st_leaves <- st_labels[!(st_labels[,8] %in% st_labels[,1]),1]

  #Identify migrations with child and parent inside subtree
  rm_mig_pos <- (is.na(ED[,9])) & (ED[,8] %in% st_labels[,1]) & (ED[,7] %in% st_labels[,1])
  ED_rm_rows <- node_indices[ED[rm_mig_pos,1]]

  internal_coal_nodes <- st_labels[!(st_labels[,1] %in% c(st_root, st_leaves)),1]

  if (is.na(st_labels[1,7])){ #If st_root is global root and not fully determined
    if (all(st_labels[1, 8:9] %in% st_labels[,1])){
      internal_coal_nodes <- append(internal_coal_nodes, st_root)
    }
  }

  int_coal_node_row <- node_indices[internal_coal_nodes]

  #### Sample new deme at internal coalescent nodes
  n_deme <- nrow(bit_rates)
  ED[int_coal_node_row, 5] <- sample(1:n_deme, length(internal_coal_nodes), TRUE)

  ### Complete migration history
  max_label <- ED[nrow(ED), 1] #max(ED[,1]) is label of final row
  log_like <- 0 #Proposal log likelihood

  for (row_id in 2 : nrow(st_labels)){ #Row 1 = st_root
    node_row <- node_indices[st_labels[row_id, 1]]
    node_deme <- ED[node_row, 5]

    parent_row <- node_indices[st_labels[row_id, 7]]
    parent_deme <- ED[parent_row, 5]

    edge_length <- ED[node_row, 6] - ED[parent_row, 6]

    trans_mat <- Re(eigen_vecs %*% diag(exp(eigen_vals * edge_length))  %*% inverse_vecs)

    mig_path <- sample_path(node_deme, parent_deme,
                            0, edge_length,
                            Q = bit_rates, P = trans_mat)

    n_mig <- nrow(mig_path) - 2
    which_child <- which(ED[parent_row, 8:9] == ED[node_row, 1])

    #log(bit_rates[i,i]) = NaN; removed by sum(..., na.rm=TRUE)
    log_like <- log_like +
      suppressWarnings(sum(log(bit_rates[mig_path[-(1 + 0:n_mig), 2] + n_deme * (mig_path[-1, 2] - 1)]), na.rm=TRUE)) + #sum(log(\lambda_{ij}))
      sum(bit_rates[mig_path[-nrow(mig_path), 2] + n_deme * (mig_path[-nrow(mig_path), 2] - 1)] * #lambda_{i+}
            (mig_path[-nrow(mig_path),1] - mig_path[-1, 1])) - # * (t_j - t_i)
      # sum(bit_rates[mig_path[-(1 + 0:n_mig), 2] + n_deme * (mig_path[-(1 + 0:n_mig), 2] - 1)] * (mig_path[-1,1] - mig_path[-(1 + 0:n_mig), 1])) - # - sum(\lambda_{i+} * (t_j - t_i))
      log(trans_mat[node_deme, parent_deme]) #Dividing by conditional probability

    if (n_mig > 0){ #Add new migrations
      new_migrations <- cbind(max_label + 1:n_mig, #New migration IDs
                              max_label + 1 + 1:n_mig, #Parent ID = next new migration ID
                              max_label - 1 + 1:n_mig, NA, #Child ID = previous new migration ID
                              mig_path[1 + 1:n_mig, 2], #New demes come from mig_path col 2 (rows 2 : nrow() - 1)
                              ED[node_row, 6] - mig_path[1 + 1:n_mig, 1],  #New times come from mig_path col 1 (rows 2 : nrow() - 1)
                              ED[parent_row, 1], #Parent coal node of current node
                              ED[node_row, 1], NA) #Child coal node of all new migrations is current node

      new_migrations[1, 3] <- ED[node_row, 1] #Correct first migration child
      new_migrations[n_mig, 2] <- ED[parent_row, 1] #Correct last migration parent

      ED[node_row, 2] <- max_label + 1 #Correct parent of current node
      ED[parent_row, 2 + which_child] <- max_label + n_mig #Correct child of parent_node

      ED <- rbind(ED, new_migrations)
      max_label <- max_label + n_mig
    } else { #Correct current and parent nodes to remove migrations
      ED[node_row, 2] <- ED[parent_row, 1] #Parent of current row is parent_row
      ED[parent_row, 2 + which_child] <- ED[node_row, 1] #Child of parent_row is current_row
    }

  }
  ### Remove old migration events
  if (length(ED_rm_rows) > 0){
    ED <- ED[-ED_rm_rows, ]
  }

  return(list(proposal=structure(ED, class='ED'), prop_prob = log_like))
}


#' MultiTypeTree Probability Density
#'
#' Computes the transition kernel for a MultiTypeTree NodeRetype operator for a full structured phylogeny
#'
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param bit_mig_mat Matrix of initial backwards-in-time migration rates for the MCMC chain
#' @param node_indices Vector of row indices corresponding to which row of ED corresponds to each node label
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to coalescent nodes in the selected subtree only (local likelihood only)
#'
#' @return List consisting of the log.likelihood and likelihood of the structured genealogy
#'
#' @export

MTT_transition_kernel <- function(ED, bit_rates, node_indices = NodeIndicesC(ED), eigen_vals = eigen(bit_rates)$values, eigen_vecs = eigen(bit_rates)$vectors, inverse_vecs = solve(eigen_vecs)){
  log_like <- 0
  for (row_id in 1 : nrow(ED)){
    if (!is.na(ED[row_id, 2])){
      parent_row <- node_indices[ED[row_id, 2]]
      edge_length <- ED[row_id, 6] - ED[parent_row, 6]

      current_deme <- ED[row_id, 5]
      parent_deme <- ED[parent_row, 5]

      log_like <- log_like + edge_length * bit_rates[current_deme, current_deme] #bit_rates[i,i] = rowSums(bit_mig_mat)[i]

      if (parent_deme != current_deme){
        log_like <- log_like + log(bit_rates[current_deme, parent_deme])
      }

      if (is.na(ED[row_id,3]) | !is.na(ED[row_id,4])){
        parent_row <- node_indices[ED[row_id, 7]]
        edge_length <- ED[row_id, 6] - ED[parent_row, 6]
        trans_prob <- Re(sum(exp(eigen_vals * edge_length) * eigen_vecs[current_deme,] * inverse_vecs[, ED[parent_row, 5]]))
        log_like <- log_like - log(trans_prob)
      }
    }
  }
  log_like <- unname(log_like)
  return(list(log.likelihood = log_like, likelihood = exp(log_like)))
}
