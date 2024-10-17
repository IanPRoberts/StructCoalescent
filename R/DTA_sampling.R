#' DTA Sampling
#'
#' Samples N migration histories under the DTA model deme fixed at all leaves using belief propagation
#'
#'
#' @param strphylo Dated structured phylogenetic tree
#' @param fit_mig_mat Forwards-in-time migration matrix for the phylogeny
#' @param N Number of migration histories to sample
#' @param parallel Logical, whether or not to generate proposals in parallel
#' @param mc.cores (optional) If parallel is TRUE, number of cores to parallelise over
#'
#' @return List of proposed DTA migration histories with the same leaves as input ED
#'
#' @export

DTA_sampling <- function(strphylo, fit_mig_mat, N = 1, parallel = FALSE, mc.cores = NA){
  if (N > 1){
    if (parallel){
      if (is.na(mc.cores)){
        options(mc.cores = parallel::detectCores() /2)
      }
      options(mc.cores = mc.cores)
    }
  }

  fit_rates <- fit_mig_mat
  diag(fit_rates) <- 0
  diag(fit_rates) <- - rowSums(fit_rates)

  ED <- as.ED(strphylo)
  topology <- ED[is.na(ED[,3]) | !is.na(ED[,4]),]
  topology[,2:4] <- topology[,7:9]
  node_indices <- NodeIndicesC(topology)
  n <- (dim(topology)[1] + 1)/2
  n_deme <- nrow(fit_rates)

  ### Transition matrices
  eigen_decomp <- eigen(fit_rates)
  V <- eigen_decomp$vectors
  V_inv <- solve(V)
  edge_lengths <- topology[,6] - topology[node_indices[topology[,2]], 6] #edge_lengths[i] = time between node at row i and its parent; edge_lengths[root] = NA

  trans_mats <- array(sapply(edge_lengths, function(x) V %*% diag(exp(eigen_decomp$values * x)) %*% V_inv),
                      dim = c(n_deme, n_deme, 2*n-1))

  ### Backwards messages
  messages <- array(0, dim = c(2*n-1, 2*n-1, n_deme))
  node_order <- order(topology[,6], decreasing = TRUE) #Ordering of coalescent nodes from lowest to root. Only need do once for entire MCMC!

  for (i in 1 : (2*n-2)){ #Final entry is root node - no messages here
    node_row <- node_order[i]
    node_parent <- topology[node_row, 2]
    node_parent_row <- node_indices[node_parent]

    if (is.na(topology[node_row, 3])){ #Leaf node
      leaf_deme <- topology[node_row, 5]
      messages[node_row, node_parent_row, ] <- trans_mats[ ,leaf_deme, node_row] #Fully determined leaf nodes
    } else{ #Coalescent node
      node_children <- topology[node_row, 3:4]
      node_children_rows <- node_indices[node_children]
      messages[node_row, node_parent_row, ] <- trans_mats[,, node_row] %*% apply(messages[node_children_rows, node_row,], 2, prod)
      messages[node_row, node_parent_row, ] <- messages[node_row, node_parent_row, ] / sum(messages[node_row, node_parent_row, ]) #Normalise messages to prevent numerical underflow
    }
  }

  proposal_func <- function(x){
    proposal <- topology #New proposal
    max_label <- max(proposal[,1])

    for (i in (2*n-1):1){
      node_row <- node_order[i]

      if (!is.na(proposal[node_row, 3])){ #Coalescent node
        node_children <- proposal[node_row, 3:4]
        node_children_rows <- node_indices[node_children]

        node_dist <- apply(messages[node_children_rows, node_row,], 2, prod) #Product of root-wards child messages

        if (!is.na(proposal[node_row, 2])){ #Non-root coalescent node
          parent_row <- node_indices[proposal[node_row,2]]
          parent_deme <- proposal[parent_row, 5]
          node_dist <- node_dist * trans_mats[parent_deme,,node_row] #Multiply by leaf-wards parent message
        }

        proposal[node_row, 5] <- sample(1:n_deme, 1, prob = node_dist) #Sample node distribution with probability proportional to node_dist
      }

      if (!is.na(proposal[node_row, 2])){ #Non-root nodes
        parent_row <- node_indices[proposal[node_row,2]] #Recalculate parent_row in case of leaf node
        mig_path <- ECctmc::sample_path(proposal[parent_row, 5], proposal[node_row, 5], 0, edge_lengths[node_row], Q = fit_rates)
        n_mig <- dim(mig_path)[1] - 2

        if (n_mig > 0){
          parent_node <- proposal[parent_row, 1]
          which_child <- which(proposal[parent_row, 3:4] == proposal[node_row])
          proposal[parent_row, 2 + which_child] <- max_label + 1

          for (k in 1 : n_mig){
            proposal <- rbind(proposal,
                              c(max_label + 1, parent_node, max_label + 2, NA, mig_path[k, 2], proposal[parent_row, 6] + mig_path[k+1, 1]))
            max_label <- max_label + 1
            parent_node <- max_label
          }
          proposal[dim(proposal)[1], 3] <- proposal[node_row, 1]
          proposal[node_row, 2] <- max_label
        }
      }
    }
    return(proposal)
  }

  if (parallel){
    proposals <- parallel::mclapply(1:N, proposal_func)
  } else {
    proposals <- lapply(1:N, proposal_func)
  }
  return(proposals)
}


# EED_local_DTA_eigen <- function(EED, fit_rates, node_indices, eigen_vals = NULL, eigen_vecs = NULL, inverse_vecs = NULL, selected_node = NULL){
#   if (is.null(selected_node)) selected_node <- sample(EED[!is.na(EED[,4]),1], 1)
#
#   if (is.null(eigen_vals) || is.null(eigen_vecs)){
#     eigen_decomp <- eigen(fit_rates)
#     eigen_vals <- eigen_decomp$values
#     eigen_vecs <- eigen_decomp$vectors
#     inverse_vecs <- solve(eigen_vecs)
#   }
#
#   if (is.null(inverse_vecs)) inverse_vecs <- solve(eigen_vecs)
#
#   selected_row <- node_indices[selected_node]
#
#   child_nodes <- EED[selected_row, 8:9]
#   parent_nodes <- rep(selected_node, 2)
#   n_edge <- 2
#
#   if (!is.na(EED[selected_row, 7])){ #If selected_node is not the root
#     child_nodes <- append(child_nodes, selected_node)
#     parent_nodes <- append(parent_nodes, EED[selected_row, 7])
#     n_edge <- 3
#   }
#
#   child_rows <- node_indices[child_nodes]
#   parent_rows <- node_indices[parent_nodes]
#   edge_lengths <- EED[child_rows, 6] - EED[parent_rows, 6]
#
#   trans_mats <- lapply(edge_lengths, function(edge_length){
#     eigen_vecs %*% diag(exp(edge_length * eigen_vals)) %*% inverse_vecs
#   })
#
#   n_deme <- nrow(fit_rates)
#
#   node_dist <- rep(1, n_deme)
#   for (x in 1 : 2){
#     node_dist <- node_dist * trans_mats[[x]][, EED[child_rows[x], 5]]
#   }
#
#   if (n_edge == 3){
#     node_dist <- node_dist * trans_mats[[3]][EED[parent_rows[3], 5], ]
#   }
#   node_dist <- node_dist / sum(node_dist)
#
#   ################################################################################
#   prop <- EED
#   prop[selected_row, 5] <- sample.int(n_deme, 1, prob = node_dist) #Update deme at selected node
#
#   # Rows to remove from proposal
#   # Any row which has exactly one of child_nodes
#   # (only node with exactly 2 is selected_node)
#   rm_rows <- which(xor(prop[,8] %in% child_nodes,
#                        prop[,9] %in% child_nodes))
#
#   if (!is.na(prop[selected_row, 7])){ #Remove parent(selected_node) from rm_rows
#     rm_rows <- rm_rows[prop[rm_rows, 1] != prop[selected_row, 7]]
#   }
#
#   child_demes <- prop[child_rows, 5]
#   child_times <- prop[child_rows, 6]
#   parent_demes <- prop[parent_rows, 5]
#   parent_times <- prop[parent_rows, 6]
#
#   max_label <- max(prop[,1])
#
#   for (x in 1 : n_edge){
#     #Provide eigendecomposition to sample_path function
#     mig_path <- ECctmc::sample_path(a = parent_demes[x],
#                                     b = child_demes[x],
#                                     t0 = parent_times[x],
#                                     t1 = child_times[x],
#                                     Q = fit_rates,
#                                     eigen_vals = eigen_vals,
#                                     eigen_vecs = eigen_vecs,
#                                     inverse_vecs = inverse_vecs)
#     n_mig <- nrow(mig_path) - 2 #Number of migration events to add
#
#     if (n_mig == 0){
#       # If no migrations between parent and child, still (might) need to update parent/child
#       prop[child_rows[x], 2] <- parent_nodes[x]
#       which_child <- which(prop[parent_rows[x], 8:9] == child_nodes[x])
#       prop[parent_rows[x], 2 + which_child] <- child_nodes[x]
#     } else {
#       # Add migrations if n_mig > 0
#       add_rows <- cbind(
#         max_label + 1:n_mig, #Node ID
#         max_label + 1:n_mig - 1, #Parent
#         max_label + 1:n_mig + 1, #Child 1
#         NA, #Child 2
#         mig_path[1:n_mig,2], # Deme
#         mig_path[1 + 1:n_mig,1], #Node Age
#         parent_nodes[x], #Parent coal
#         child_nodes[x], #Child coal 1
#         NA #Child coal 2
#       )
#
#       # Update first and final to have child/parent_nodes[x] in correct place
#       add_rows[1, 2] <- parent_nodes[x]
#       add_rows[n_mig, 3] <- child_nodes[x]
#
#       prop <- rbind(prop, add_rows)
#
#       # Update child/parent_nodes[x]
#       prop[child_rows[x], 2] <- max_label + n_mig
#       which_child <- which(prop[parent_rows[x], 8:9] == child_rows[x])
#       prop[parent_rows[x], 2 + which_child] <- max_label + 1
#     }
#     # Update max_label
#     max_label <- max_label + n_mig
#   }
#
#   if (length(rm_rows) > 0) prop <- prop[-rm_rows,]
#   return(list(proposal = prop, updated_node = selected_node, node_dist = node_dist))
# }



local_DTA_subtree_proposal <- function(EED, st_labels, fit_rates, node_indices = NodeIndicesC(EED), eigen_decomp = eigen(fit_rates), inverse_vecs = solve(eigen_decomp$vectors)){
  n_deme <- nrow(fit_rates)

  #Identify migrations with child and parent inside subtree
  st_rm_rows <- is.na(st_labels[,4]) & (st_labels[,2] %in% st_labels[,1]) & (st_labels[,3] %in% st_labels[,1])
  EED_rm_rows <- node_indices[st_labels[st_rm_rows, 1]]

  #Remove migrations from inside st_labels
  st_labels <- st_labels[!st_rm_rows,]

  ######################
  # Backwards messages #
  ######################
  #Reorder subtree nodes into decreasing node age
  st_order <- order(st_labels[,6], decreasing = TRUE)
  st_labels <- st_labels[st_order,]
  st_rows <- node_indices[st_labels[,1]]


  n_nodes <- nrow(st_labels)
  st_root <- st_labels[n_nodes,1] #Label of st_root -> node with least node age
  st_root_row <- node_indices[st_root]

  #Parent coal nodes within subtree
  st_parent_labels <- EED[st_rows, 7]
  st_parent_labels[st_parent_labels == st_parent_labels[n_nodes]] <- st_root
  st_parent_labels[n_nodes] <- st_root #In case st_root is global root
  st_parent_rows <- node_indices[st_parent_labels]
  st_edge_lengths <- EED[st_rows, 6] - EED[st_parent_rows, 6]
  st_edge_lengths[is.na(st_edge_lengths)] <- 0


  #Child coal nodes within subtree
  st_children <- matrix(NA, n_nodes, 2) #unname(EED[st_rows, 8:9])
  st_child_rows <- matrix(NA, n_nodes, 2) #matrix(node_indices[st_children], ncol = 2)

  for (st_id in 1 : n_nodes){
    for (child_id in 1 : 2){
      child_label <- EED[st_rows[st_id], 7 + child_id]
      child_row <- node_indices[child_label]

      if (!is.na(child_row)){
        while (!(child_label %in% st_labels[,1])){
          child_label <- EED[child_row, 2]
          child_row <- node_indices[child_label]
        }
      }
      st_children[st_id, child_id] <- child_label
      st_child_rows[st_id, child_id] <- child_row
    }
  }

  st_child_ids <- matrix(NA, n_nodes, 2)
  st_parent_ids <- numeric(n_nodes)
  for (st_id in 1 : n_nodes){
    st_parent_ids[st_parent_labels == st_labels[st_id]] <- st_id
    st_child_ids[st_children == st_labels[st_id]] <- st_id
  }

  #Transition matrices
  trans_mats <- array(0, c(n_deme, n_deme, n_nodes))
  for (st_id in 1:(n_nodes - 1)){ #Don't need to exponentiate to identity matrix for edge length 0 above root
    trans_mats[,,st_id] <- Re(eigen_decomp$vectors %*% diag(exp(st_edge_lengths[st_id] * eigen_decomp$values)) %*% inverse_vecs)
  }

  messages <- array(NA, dim = c(n_nodes, n_nodes, n_deme))

  for (st_id in 1 : (n_nodes - 1)){ #Skip subtree root (st_root corresponds to st_id == n_nodes)
    node_row <- st_rows[st_id]

    if ((is.na(st_labels[st_id, 3])) || st_children[st_id, 1] == EED[node_row, 1]){
      # Subtree leaf node contributes column of transition matrix ending at current deme
      current_deme <- EED[node_row, 5]
      messages[st_id, st_parent_ids[st_id], ] <- trans_mats[, current_deme, st_id]
    } else {
      # Other subtree nodes contribute transition matrix %*% incoming messages
      # messages[st_id, st_parent_ids[st_id], ] <- trans_mats[,, st_id] %*% apply(messages[st_child_ids[st_id,], st_id,], 2, prod)
      messages[st_id, st_parent_ids[st_id], ] <- trans_mats[,, st_id] %*% (messages[st_child_ids[st_id,1], st_id ,] * messages[st_child_ids[st_id,2], st_id, ])
      messages[st_id, st_parent_ids[st_id], ] <- messages[st_id, st_parent_ids[st_id], ] / sum(messages[st_id, st_parent_ids[st_id], ])
    }
  }

  #####################
  # Forwards sampling #
  #####################
  #When st_root is global root (a.s. the only situation where st_root may be coalescent node)
  #Need to sample deme at st_root as well as internal node!

  if (is.na(EED[st_rows[n_nodes], 2])){ #st_root == global_root
    #Deme distribution just product of incoming messages from below
    node_dist <- messages[st_child_ids[n_nodes, 1], n_nodes,] * messages[st_child_ids[n_nodes, 2], n_nodes,]
    EED[st_root_row, 5] <- sample.int(n_deme, 1, prob = node_dist)
  }

  for (st_id in (n_nodes - 1):1){ #Loop over non-subtree-leaf coalescent nodes
    node_row <- st_rows[st_id]
    if ((!is.na(st_child_rows[st_id, 1])) && (st_child_rows[st_id, 1] != node_row)){
      parent_deme <- EED[st_parent_rows[st_id], 5]

      #Deme distribution is product of transition from determined parent with product of backwards messages from node's children
      node_dist <- trans_mats[parent_deme,,st_id] * messages[st_child_ids[st_id,1], st_id,] * messages[st_child_ids[st_id,2], st_id,]
      EED[node_row, 5] <- sample.int(n_deme, 1, prob = node_dist)
    }
  }

  ##################
  # Subtree update #
  ##################
  max_label <- max(EED[,1])
  log_like <- 0 #Proposal log likelihood

  #Fill in parent edge of each node in st_labels except st_root
  for (st_id in 1 : (n_nodes - 1)){ #Skip subtree root (a.s. last in node_order)
    node_row <- st_rows[st_id]
    node_label <- st_labels[st_id]

    parent_deme <- EED[st_parent_rows[st_id], 5]
    parent_time <- EED[st_parent_rows[st_id], 6]

    # mig_path <- ECctmc::sample_path_unif(parent_deme, EED[node_row, 5], 0, st_edge_lengths[st_id], Q = fit_rates)
    mig_path <- ECctmc::sample_path_unif3(parent_deme, EED[node_row, 5], 0, st_edge_lengths[st_id],
                                          Q = fit_rates, P = trans_mats[,,st_id])


    ###### ECctmc::sample_path_unif returns t(mig_path) if 0 migrations added
    ###### BUG IN SOURCE CODE!!!
    if (nrow(mig_path) == ncol(mig_path)){
      if (mig_path[1,2] == st_edge_lengths[st_id]){
        mig_path <- t(mig_path)
      }
    }

    n_mig <- nrow(mig_path) - 2

    which_child <- which(st_children[st_parent_ids[st_id], ] == node_label)

    if (n_mig > 0){ #Add new migrations
      parent_node <- st_parent_labels[st_id]

      EED[st_parent_rows[st_id], 2 + which_child] <- max_label + 1 #Update child of node_parent
      for (k in 1 : n_mig){
        #log-probability of holding time until next migration i -> j
        log_like <- log_like +
          log(fit_rates[mig_path[k, 2], mig_path[k + 1, 2]]) + #log(f_ij)
          fit_rates[mig_path[k, 2], mig_path[k, 2]] * (mig_path[k+1, 1] - mig_path[k,1]) # - f_{i+} * (t_j - t_i)

        max_label <- max_label + 1
        EED <- rbind(EED,
                     c(max_label, #New migration ID
                       parent_node,
                       max_label + 1, NA, #Next migration ID
                       mig_path[k, 2], #New migration deme
                       parent_time + mig_path[k+1, 1], #Migration time increases leaf-wards
                       EED[node_row, 7], #Maintain same coal node parent as node_row
                       EED[st_parent_rows[st_id], 7 + which_child], NA))
        parent_node <- max_label
        node_indices[max_label] <- nrow(EED)
      }
      EED[nrow(EED), 3] <- node_label #Update child of final migration added
      EED[node_row, 2] <- max_label #Update parent of current node
    } else{ #No migrations added
      EED[node_row, 2] <- st_parent_labels[st_id] #Update parent of current node
      EED[st_parent_rows[st_id], 2 + which_child] <- node_label #Update child of parent node
    }

    #log-probability of no further migrations between last two events on mig_path
    log_like <- log_like + fit_rates[mig_path[n_mig + 1, 2], mig_path[n_mig + 2, 2]] * (mig_path[n_mig+2, 1] - mig_path[n_mig + 1, 1])
  }

  #############################
  # Remove virtual migrations #
  ############################
  # If st_root is a migration, it is a.s. self-migration
  # Else st_root is the global root and no change needs to be made
  if (is.na(EED[st_root_row, 4])){
    parent_node <- EED[st_root_row, 2]
    parent_row <- node_indices[parent_node]
    which_child <- which(EED[parent_row, 8:9] == EED[st_root_row, 8])
    EED[parent_row, 2 + which_child] <- EED[st_root_row, 3] #Child of parent is now child of st_root

    child_node <- EED[st_root_row, 3]
    child_row <- node_indices[child_node]
    EED[child_row, 2] <- parent_node

    EED_rm_rows <- append(EED_rm_rows, st_root_row)
  }

  for (st_id in 1 : (n_nodes - 1)){
    current_row <- st_rows[st_id]

    if ((is.na(EED[current_row, 4])) & (!is.na(EED[current_row, 3]))){ #If current_row is migration then a.s. self-migration
      child_node <- EED[current_row, 3]
      child_row <- node_indices[child_node]
      EED[child_row, 2] <- EED[current_row, 2] #Parent of child_node is parent(current_st_leaf)

      parent_node <- EED[current_row, 2]
      parent_row <- node_indices[parent_node]
      which_child <- which(EED[parent_row, 3:4] == st_labels[st_id, 1])
      EED[parent_row, 2 + which_child] <- child_node

      EED_rm_rows <- append(EED_rm_rows, current_row)
    }
  }

  EED <- EED[-EED_rm_rows,]

  return(list(proposal=EED, prop_prob = log_like))
}


#' Subtree Sampling
#'
#' Samples a subtree up to patristic distance st_radius from a subtree centre
#'
#'
#' @param ED ED representation of structured phylogenetic tree
#' @param st_radius Subtree radius
#' @param NI Vector of node indices. Element i gives the row of label i
#' @param st_child Child node of the subtree centre
#' @param st_centre_loc Value between 0 and 1 giving position between st_child and ED[st_child,2]
#'
#' @return List of proposed DTA migration histories with the same leaves as input ED

st_centre_dist <- function(ED, st_radius, NI, st_child = NA, st_centre_loc = runif(1)){
  root_node <- ED[is.na(ED[,2]), 1]
  root_row <- NI[root_node]

  edge_lengths <- ED[,6] - ED[NI[ED[,2]], 6]
  edge_lengths[root_row] <- 0


  if (is.na(st_child)) st_child <- sample(ED[,1], 1, prob = edge_lengths)

  st_child_row <- NI[st_child]
  st_root <- st_child
  st_root_row <- st_child_row

  st_centre_age <- ED[st_child_row, 6] - edge_lengths[st_child_row] * st_centre_loc

  st_root_age <- max(ED[root_row, 6], st_centre_age - st_radius)

  #st_labels contains node_ID, parent, children and distance from st_centre
  st_labels <- matrix(NA, nrow = 0, ncol = 2,
                      dimnames = list(NULL, c("Node_ID", "Node_dist")))

  while (ED[st_root_row, 6] > st_root_age){
    st_labels <- rbind(st_labels,
                       ED[st_root_row, c(1, 6)])
    st_root <- ED[st_root_row, 2]
    st_root_row <- NI[st_root]
  }

  if (ED[st_root_row, 6] < st_root_age){
    #Add virtual migration as st_root if needed
    max_label <- max(ED[,1]) + 1
    which_child <- which(ED[st_root_row, 3:4] == st_labels[nrow(st_labels), 1])
    st_root_child <- ED[st_root_row, 2 + which_child]

    if (is.na(ED[st_root_row, 4])){ #If st_root is a migration, parent coal is parent coal of st_root
      parent_coal <- ED[st_root_row, 7]
    } else { #Else parent coal is st_root
      parent_coal <- st_root
    }

    ED <- rbind(ED,
                c(max_label, #Label
                  st_root, #Parent
                  ED[st_root_row, 2 + which_child], #Child 1
                  NA, #Child 2
                  ED[NI[st_root_child], 5], #Deme
                  st_root_age, #Node age
                  parent_coal, #Parent coal
                  ED[st_root_row, 7 + which_child], #Child coal 1
                  NA #Child coal 2
                ))

    ED[NI[st_root_child], 2] <- max_label
    ED[st_root_row, 2 + which_child] <- max_label
    NI[max_label] <- nrow(ED)

    st_labels <- rbind(st_labels,
                       ED[nrow(ED), c(1, 6)])
  } else { #No virtual migration required as st_root is global root
    max_label <- max(ED[,1])
    st_labels <- rbind(st_labels,
                       ED[st_root_row, c(1, 6)])
  }

  st_labels[, 2] <- st_centre_age - st_labels[, 2] #abs(st_centre_age - st_labels[, 2])

  current_parents <- st_labels[,1]

  if (st_labels[1,2] < -st_radius){
    current_parents <- current_parents[-1]
  }

  while (length(current_parents) > 0){
    new_parents <- numeric(0)
    for (node in current_parents){
      node_row <- NI[node]
      which_child <- which(!(na.omit(ED[node_row, 3:4]) %in% st_labels[,1]))
      node_dist <- st_labels[which(st_labels[,1] == node), 2]

      if (length(which_child) > 0){
        for (child_id in which_child){
          child_row <- NI[ED[node_row, 2 + child_id]]

          #### Need to remove new_parents below st_leaf_age and global leaves
          # child_dist <- node_dist + edge_lengths[child_row]
          child_dist <- node_dist - edge_lengths[child_row]
          st_labels <- rbind(st_labels,
                             c(ED[child_row, 1], child_dist))

          if ((abs(child_dist) < st_radius) & (!is.na(ED[child_row, 3]))){
            new_parents <- append(new_parents, ED[child_row, 1])
          }
        }
      }
    }
    current_parents <- new_parents
  }

  for (row_id in 1 : nrow(st_labels)){
    # if (st_labels[row_id, 2] > st_radius){
    if (st_labels[row_id, 2] < - st_radius){
      #If node is greater than distance st_radius from st_centre add virtual migration
      max_label <- max_label + 1

      node_row <- NI[st_labels[row_id, 1]]

      if ((is.na(ED[node_row, 4])) & (!is.na(ED[node_row, 3]))){
        #Current node is a migration -> child coal is child coal of current node
        child_coal <- ED[node_row, 8]
      } else {
        #Current node is a coalescent or leaf -> child coal is current node
        child_coal <- st_labels[row_id, 1]
      }

      ED <- rbind(ED,
                  c(max_label, #Label
                    ED[node_row, 2], #Parent
                    ED[node_row, 1], #Child 1
                    NA, #Child 2
                    ED[node_row, 5], #Deme
                    ED[node_row, 6] + st_radius - abs(st_labels[row_id, 2]), #Node age
                    ED[node_row, 7], #Parent coal
                    child_coal, #Child coal 1
                    NA #Child coal 2
                  ))

      parent_row <- NI[ED[node_row, 2]]
      ED[node_row, 2] <- max_label
      which_child <- which(ED[parent_row, 3:4] == st_labels[row_id, 1])
      ED[parent_row, 2 + which_child] <- max_label

      NI[max_label] <- nrow(ED)

      st_labels[row_id,] <- c(ED[nrow(ED), 1], - st_radius)
    }
  }

  return(list(ED = ED, st_labels = ED[NI[st_labels[,1]],]))
}
