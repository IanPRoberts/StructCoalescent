#' Migration Birth MCMC Move
#'
#' Performs a Migration birth move (Ewing et al. 2004). Adds a migration
#' node between a randomly selected ancestral node and its parent, allocating a
#' deme for the new edge consistent with the surrounding edges.
#'
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n_deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration birth move

migration_birth <- function(ED, n_deme, fix_leaf_deme = TRUE, node_indices){
  edge_lengths <- ED[,6] - ED[node_indices[ED[,2]], 6]
  edge_lengths[is.na(edge_lengths)] <- 0
  child_row <- node_indices[sample(ED[,1], 1, prob=edge_lengths)] #Sample child_row with probability proportional to branch length above

  subtree_rows <- active_rows <- child_row

  #Subtree containing <new_node, child_node> terminating in migration or leaf nodes
  while (length(active_rows) > 0){
    active_rows <- active_rows[!is.na(ED[active_rows, 4])] #Coalescent nodes only
    active_rows <- node_indices[ED[active_rows, 3:4]]
    subtree_rows <- c(subtree_rows,
                      active_rows)
  }

  if (fix_leaf_deme == TRUE & anyNA(ED[subtree_rows, 3])){
    #Early rejection if leaf in subtree and leaf demes are fixed
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  subtree_leaf_rows <- subtree_rows[is.na(ED[subtree_rows, 4])]
  current_deme <- ED[child_row, 5]

  if (n_deme == 2){
    new_deme <- 3 - current_deme
  } else {
    new_deme <- sample((1:n_deme)[-current_deme], 1)
  }

  subtree_leaf_children_rows <- node_indices[na.omit(ED[subtree_leaf_rows, 3])]

  if (any(ED[subtree_leaf_children_rows, 5] == new_deme)){
    #Early rejection if leaf in subtree matches new deme
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  ED[subtree_rows, 5] <- new_deme
  node_ID <- max(ED[,1]) + 1

  #Update child of parent of child_row
  parent_node <- ED[child_row, 2]
  child_node <- ED[child_row, 1]
  parent_row <- node_indices[parent_node]
  ED[parent_row, 2 + which(ED[parent_row, 3:4] == child_node)] <- node_ID

  ED[child_row, 2] <- node_ID #Update parent of child_row


  migration_age <- ED[child_row, 6] - runif(1) * edge_lengths[child_row] #Age of new migration
  if (is.na(ED[child_row, 4])){ #Child_node is migration or leaf
    ED <- rbind(ED,
                c(node_ID, parent_node, child_node, NA, current_deme, migration_age, ED[child_row, 7:9]))
  } else { #Child_node is coalescence event
    ED <- rbind(ED,
                c(node_ID, parent_node, child_node, NA, current_deme, migration_age, ED[child_row, 7], child_node, NA))
  }

  tree_length <- sum(edge_lengths, na.rm=TRUE)
  n_mig <- sum(is.na(ED[,4]) & !is.na(ED[,3]))
  prop_ratio <- (n_deme - 1) * tree_length / (n_mig + 1)
  log_prop_ratio <- log(n_deme - 1) + log(tree_length) - log(n_mig + 1)

  return(list(ED = structure(ED, class='ED'), prop.ratio=prop_ratio, log.prop.ratio=log_prop_ratio, node_indices=NodeIndicesC(ED)))
}

#' Migration Death MCMC Move
#'
#' Performs a Migration death move (Ewing et al. 2004). Selects an ancestral node
#' uniformly at random and removes the parent of the selected node, allocating
#' demes as necessary for modified edges.
#'
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n_deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration death move

migration_death <- function(ED, n_deme, fix_leaf_deme = TRUE, node_indices){
  migration_nodes <- ED[is.na(ED[,4]) & !is.na(ED[,3]), 1]

  if (length(migration_nodes) == 0){
    #Reject, no migration events in ED to remove
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  if (length(migration_nodes) == 1){
    selected_node <- migration_nodes
  } else {
    selected_node <- sample(migration_nodes, 1)
  }

  selected_row <- node_indices[selected_node]
  child_node <- ED[selected_row, 3]
  child_row <- node_indices[child_node]

  subtree_rows <- active_rows <- child_row

  #Subtree containing <selected_node, child_node> terminating in migration or leaf nodes
  while (length(active_rows) > 0){
    active_rows <- active_rows[!is.na(ED[active_rows, 4])] #Coalescent nodes only
    active_rows <- node_indices[ED[active_rows, 3:4]]
    subtree_rows <- c(subtree_rows,
                      active_rows)
  }

  if (fix_leaf_deme == TRUE & anyNA(ED[subtree_rows, 3])){
    #Early rejection if leaf in subtree and leaf demes are fixed
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  subtree_leaf_rows <- subtree_rows[is.na(ED[subtree_rows, 4])]
  subtree_leaf_children_rows <- node_indices[na.omit(ED[subtree_leaf_rows, 3])]

  proposal_deme <- ED[selected_row, 5] #Deme to propagate through subtree

  if (proposal_deme %in% ED[subtree_leaf_children_rows, 5]){
    #Reject, deme invalid
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  #Update deme of subtree
  ED[subtree_rows, 5] <- proposal_deme

  #Update parent of child_row
  parent_node <- ED[selected_row, 2]
  ED[child_row, 2] <- parent_node

  #Update child of parent_node
  parent_row <- node_indices[parent_node]
  ED[parent_row, 2 + which(ED[parent_row, 3:4] == ED[selected_row, 1])] <- child_node

  #Update node_indices
  node_indices[selected_node] <- 0
  node_indices[node_indices > selected_row] <- node_indices[node_indices > selected_row] - 1

  ED <- ED[-selected_row,]
  tree_length <- sum(ED[,6] - ED[node_indices[ED[,2]], 6], na.rm=TRUE)

  prop_ratio <- length(migration_nodes)/ ((n_deme - 1) * tree_length)
  log_prop_ratio <- log(length(migration_nodes)) - log(n_deme - 1) - log(tree_length)
  return(list(ED = structure(ED, class='ED'), prop.ratio=prop_ratio, log.prop.ratio=log_prop_ratio, node_indices=node_indices))
}

#' Migration Pair Birth MCMC Move
#'
#' Performs a Migration pair birth move (Ewing et al. 2004). Adds two migration
#' nodes on an edge selected uniformly at random from a structured coalescent
#' process, allocating a deme for the added edge such that a migration event
#' does not target its origin deme.
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n_deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration pair birth move

migration_pair_birth <- function(ED, n_deme, node_indices){
  selected_node <- sample(ED[!is.na(ED[,2]), 1], 1)
  selected_row <- node_indices[selected_node]

  parent_node <- ED[selected_row, 2]
  parent_row <- node_indices[parent_node]

  node_IDs <- max(ED[,1]) + 1:2
  node_ages <- sort(runif(2, ED[parent_row, 6], ED[selected_row, 6]), decreasing=TRUE)

  if (n_deme == 2){
    new_deme <- 3 - ED[selected_row, 5]
  } else {
    new_deme <- sample((1:n_deme)[-ED[selected_row, 5]], 1)
  }

  #Parent of selected_node
  ED[selected_row, 2] <- node_IDs[1]

  #Child of parent_node
  which_child <- which(ED[parent_row, 3:4] == selected_node)
  ED[parent_row, 2 + which_child] <- node_IDs[2]

  #New migration events
  ED <- rbind(ED,
              c(node_IDs[1], node_IDs[2], selected_node, NA, new_deme, node_ages[1], ED[selected_row, 7], ED[parent_row, 7 + which_child], NA),
              c(node_IDs[2], parent_node, node_IDs[1], NA, ED[selected_row, 5], node_ages[2], ED[selected_row, 7], ED[parent_row, 7 + which_child], NA))

  n <- nrow(ED)
  prop_ratio <- (n_deme - 1) * (n - 1) * (ED[selected_row, 6] - ED[parent_row, 6])^2 / (2 * (n + 1))
  log_prop_ratio <- log(n_deme - 1) + log(n - 1) + 2 * log(ED[selected_row, 6] - ED[parent_row, 6]) - log(2) - log(n + 1)

  if (node_IDs[2] > length(node_indices)){
    new_node_indices <- numeric(node_IDs[2])
    new_node_indices[1:length(node_indices)] <- node_indices
  } else{
    new_node_indices <- node_indices
  }
  new_node_indices[node_IDs] <- n - 1:0
  return(list(ED = structure(ED, class='ED'), prop.ratio = prop_ratio, log.prop.ratio = log_prop_ratio, node_indices = new_node_indices))
}

#' Migration Pair Death MCMC Move
#'
#' Performs a Migration pair death move (Ewing et al. 2004). Deletes two
#' migration nodes if they lie between two edges in the same deme.
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n_deme Number of distinct demes in the population
#'
#' @return Updated extended data object with the proposal from the migration pair death move

migration_pair_death <- function(ED, n_deme, node_indices){
  selected_node <- sample(ED[!is.na(ED[,2]), 1], 1)
  selected_row <- node_indices[selected_node]

  if (is.na(ED[selected_row, 3]) | !is.na(ED[selected_row, 4])){
    #Reject, selected_node not a migration event
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  parent_node <- ED[selected_row, 2]
  parent_row <- node_indices[parent_node]

  if (is.na(ED[parent_row, 3]) | !is.na(ED[parent_row, 4])){
    #Reject, parent_node not a migration event
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  child_row <- node_indices[ED[selected_row, 3]]
  child_deme <- ED[child_row, 5]

  if (child_deme == ED[parent_row, 5]){
    #Update child of parent(parent_node)
    parent_parent_row <- node_indices[ED[parent_row, 2]]
    which_child <- which(ED[parent_parent_row, 3:4] == parent_node)
    ED[parent_parent_row, 2 + which_child] <- ED[child_row, 1]

    #Update parent of child_node
    ED[child_row, 2] <- ED[parent_parent_row, 1]

    n <- nrow(ED)
    prop_ratio <- 2 * (n-1) / ((n-3) * (n_deme - 1) * (ED[child_row, 6] - ED[parent_parent_row, 6])^2)
    log_prop_ratio <- log(2) + log(n-1) - log(n-3) - log(n_deme - 1) - 2 * log(ED[child_row, 6] - ED[parent_parent_row, 6])

    ED <- ED[-c(selected_row, parent_row),]
  } else {
    #Reject, exterior demes don't match
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  return(list(ED = structure(ED, class='ED'), prop.ratio=prop_ratio, log.prop.ratio=log_prop_ratio, node_indices=NodeIndicesC(ED)))
}

#' Coalescent Node Merge Proposal
#'
#' Performs a coalescent node merge move (Ewing et al. 2004). Merges two migration
#' nodes immediately below a coalescent node and place above coalescent node
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n_deme Number of distinct demes in the population
#' @param node_indices Vector giving row indices for node labels
#'
#' @return Updated extended data object with the proposal from the migration pair birth move

coalescent_merge <- function(ED, n_deme, node_indices){
  selected_node <- sample(ED[!is.na(ED[,4]), 1], 1)
  selected_row <- node_indices[selected_node]

  child_nodes <- ED[selected_row, 3:4]
  child_rows <- node_indices[child_nodes]

  if (any(!is.na(ED[child_rows, 4])) | any(is.na(ED[child_rows, 3]))){
    #Reject, at least one child is coalescent or leaf node
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  child_child_nodes <- ED[child_rows, 3]
  child_child_rows <- node_indices[child_child_nodes]
  child_demes <- ED[child_child_rows, 5]

  if (child_demes[1] != child_demes[2]){
    #Reject, child demes not identical
    return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
  }

  if (is.na(ED[selected_row, 2])){
    #Selected_node is root - no parent migration added
    prop_ratio <- 1 / ((n_deme - 1) * prod(ED[child_child_rows, 6] - ED[selected_row, 6]))
    log_prop_ratio <- - log(n_deme - 1) - sum(log(ED[child_child_rows, 6] - ED[selected_row, 6]))

    #Update selected_row
    ED[selected_row, 3:4] <- child_child_nodes
    ED[selected_row, 5] <- child_demes[1]
  } else {
    #Selected_node is not root - parent migration added
    parent_node <- ED[selected_row, 2]
    parent_row <- node_indices[parent_node]

    node_ID <- max(ED[,1]) + 1
    node_age <- runif(1, ED[parent_row, 6], ED[selected_row, 6])

    #Add new node
    ED <- rbind(ED,
                c(node_ID, parent_node, selected_node, NA, ED[selected_row, 5], node_age, ED[selected_row, 7], selected_node, NA))

    #Update selected_node
    ED[selected_row, 2:5] <- c(node_ID, child_child_nodes, child_demes[1])

    #Update child_child_nodes
    ED[child_child_rows, 2] <- selected_node

    #Update parent_node
    which_child <- which(ED[parent_row, 3:4] == selected_node)
    ED[parent_row, 2 + which_child] <- node_ID

    prop_ratio <- (ED[selected_row, 6] - ED[parent_row, 6]) / prod(ED[child_child_rows, 6] - ED[selected_row, 6])
    log_prop_ratio <- log(ED[selected_row, 6] - ED[parent_row, 6]) - sum(log(ED[child_child_rows, 6] - ED[selected_row, 6]))
  }

  #Remove child_rows
  ED <- ED[-child_rows,]
  return(list(ED = structure(ED, class='ED'), prop.ratio=prop_ratio, log.prop.ratio=log_prop_ratio, node_indices=NodeIndicesC(ED)))
}

#' Coalescent Node Split Proposal
#'
#' Performs a coalescent node split move (Ewing et al. 2004). Splits a migration
#' node immediately above a coalescent node into two migration nodes placed immediately
#' below the coalescent node
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n_deme Number of distinct demes in the population
#' @param node_indices Vector giving row indices for node labels
#'
#' @return Updated extended data object with the proposal from the migration pair birth move

coalescent_split <- function(ED, n_deme, node_indices){
  selected_node <- sample(ED[!is.na(ED[,4]), 1], 1)
  selected_row <- node_indices[selected_node]

  child_nodes <- ED[selected_row, 3:4]
  child_rows <- node_indices[child_nodes]

  node_IDs <- max(ED[,1]) + 1:2
  node_ages <- sapply(ED[child_rows, 6], runif, n=1, min=ED[selected_row, 6])

  current_deme <- ED[child_rows[1], 5]
  parent_node <- ED[selected_row, 2]

  if (is.na(parent_node)){
    #Selected_node is root - no check for parent migration
    if (n_deme == 2){
      proposal_deme <- 3 - current_deme
    } else {
      proposal_deme <- sample((1:n_deme)[-current_deme], 1)
    }

    prop_ratio <- (n_deme - 1) * prod(ED[child_rows, 6] - ED[selected_row, 6])
    log_prop_ratio <- log(n_deme - 1) + sum(log(ED[child_rows, 6] - ED[selected_row, 6]))
  } else {
    #Selected_node is not root - check for parent migration
    parent_row <- node_indices[parent_node]

    if (is.na(ED[parent_row, 4])){
      #Reject, parent_node is coalescence
      return(list(ED = structure(ED, class='ED'), prop.ratio=0, log.prop.ratio=-Inf, node_indices=node_indices))
    }

    parent_parent_node <- ED[parent_row, 2]
    parent_parent_row <- node_indices[parent_parent_node]

    proposal_deme <- ED[selected_row, 5]

    prop_ratio <- prod(ED[child_rows, 6] - ED[selected_row, 6]) / (ED[selected_row, 6] - ED[parent_parent_row, 6])
    log_prop_ratio <- sum(log(ED[child_rows, 6] - ED[selected_row, 6])) - log(ED[selected_row, 6] - ED[parent_parent_row, 6])

    #Update children of parent_parent_node
    which_child <- which(ED[parent_parent_row, 3:4] == parent_node)
    ED[parent_parent_row, 2 + which_child] <- selected_node

    #Update parent of selected_node
    ED[selected_row, 2] <- parent_parent_node
  }

  #Update selected_node
  ED[selected_row, 3:4] <- node_IDs
  ED[selected_row, 5] <- proposal_deme

  #Update child_nodes
  ED[child_rows, 2] <- node_IDs

  #Add new nodes
  ED <- rbind(ED,
              c(node_IDs[1], selected_node, child_nodes[1], NA, proposal_deme, node_ages[1], ED[child_rows[1], 7], ED[selected_row, 8], NA),
              c(node_IDs[2], selected_node, child_nodes[2], NA, proposal_deme, node_ages[2], ED[child_rows[2], 7], ED[selected_row, 8], NA))

  if (!is.na(parent_node)){
    #Remove parent_node
    ED <- ED[-parent_row,]
  }

  return(list(ED = structure(ED, class='ED'), prop.ratio=prop_ratio, log.prop.ratio=log_prop_ratio, node_indices=NodeIndicesC(ED)))
}
