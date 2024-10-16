#' MultiTypeTree Subtree
#'
#' Extracts a subtree (and associated sub-migration history) of a structured
#' genealogy consisting of the parent branch of a coalescent event and a fixed
#' number of generations of descendant branches
#'
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param st_depth Number of generations below the selected coalescent node to include in the subtree
#' @param NI Vector of row indices corresponding to which row of ED corresponds to each node label
#' @param selected_node (optional) Selected coalescent node to use as centre of subtree
#'
#' @return List consisting of ED, the structured phylogeny input, and st_labels, a reduced extended data structure consisting of the coalescent nodes in the subtree only
#'
#' @export

MTT_st_coal_node <- function(ED, st_depth = 1, NI = NodeIndicesC(ED), selected_node = NA){
  if (is.na(selected_node)){
    coal_nodes <- ED[!is.na(ED[,4]), 1]
    selected_node <- sample(coal_nodes, 1)
  }

  selected_row <- NI[selected_node]
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
    active_rows <- NI[child_coals]
  }


  return(list(ED = ED, st_labels = ED[NI[st_labels],]))
}

#' MultiTypeTree Node Retype
#'
#' Generates a new migration history on a subtree drawn using MTT_st_coal_node()
#'
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to coalescent nodes in the selected subtree only
#' @param bit_rates Transition matrix of a Markov process corresponding to a process with rates given by the backwards-in-time migration matrix
#' @param NI Vector of row indices corresponding to which row of ED corresponds to each node label
#'
#' @return  List consisting of proposal, a structured phylogeny with updates made to the selected subtree, and prop_prob, the probability of the given update
#'
#' @export

MTT_node_retype <- function(ED, st_labels, bit_rates, NI = NodeIndicesC(ED)){
  st_root <- unname(st_labels[1,1]) #st_labels[1,] always corresponds to st_root by construction
  st_leaves <- st_labels[!(st_labels[,8] %in% st_labels[,1]),1]

  #Identify migrations with child and parent inside subtree
  rm_mig_pos <- (is.na(ED[,9])) & (ED[,8] %in% st_labels[,1]) & (ED[,7] %in% st_labels[,1])
  ED_rm_rows <- NI[ED[rm_mig_pos,1]]

  internal_coal_nodes <- st_labels[!(st_labels[,1] %in% c(st_root, st_leaves)),1]

  if (is.na(st_labels[1,7])){ #If st_root is global root and not fully determined
    if (all(st_labels[1, 8:9] %in% st_labels[,1])){
      internal_coal_nodes <- append(internal_coal_nodes, st_root)
    }
  }

  int_coal_node_row <- NI[internal_coal_nodes]

  #### Sample new deme at internal coalescent nodes
  n_deme <- nrow(bit_rates)
  ED[int_coal_node_row, 5] <- sample(1:n_deme, length(internal_coal_nodes), TRUE)

  ### Complete migration history
  max_label <- max(ED[,1])
  log_like <- 0 #Proposal log likelihood

  for (row_id in 2 : nrow(st_labels)){ #Row 1 = st_root
    node_row <- NI[st_labels[row_id, 1]]
    node_deme <- ED[node_row, 5]

    parent_row <- NI[st_labels[row_id, 7]]
    parent_deme <- ED[parent_row, 5]

    edge_length <- ED[node_row, 6] - ED[parent_row, 6]

    mig_path <- ECctmc::sample_path_unif(node_deme, parent_deme,
                                         0, edge_length,
                                         bit_rates)

    ###### ECctmc::sample_path_unif returns t(mig_path) if 0 migrations added
    ###### BUG IN SOURCE CODE!!!
    if (nrow(mig_path) == ncol(mig_path)){
      if (mig_path[1,2] == edge_length){
        mig_path <- t(mig_path)
      }
    }

    n_mig <- nrow(mig_path) - 2
    which_child <- which(ED[parent_row, 8:9] == ED[node_row, 1])

    trans_mat <- expm::expm(bit_rates * edge_length)

    #log(bit_rates[i,i]) = NaN; removed by sum(..., na.rm=TRUE)
    log_like <- log_like +
      suppressWarnings(sum(log(bit_rates[mig_path[-(1 + 0:n_mig), 2] + n_deme * (mig_path[-1, 2] - 1)]), na.rm=TRUE)) + #sum(log(\lambda_{ij}))
      sum(bit_rates[mig_path[-(1 + 0:n_mig), 2] + n_deme * (mig_path[-(1 + 0:n_mig), 2] - 1)] * (mig_path[-1,1] - mig_path[-(1 + 0:n_mig), 1])) - # - sum(\lambda_{i+} * (t_j - t_i))
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

  return(list(proposal=ED, prop_prob = log_like))
}

#' @rdname MTT_node_retype
#' @export

MTT_node_retype_eigen <- function(ED, st_labels, bit_rates, NI = NodeIndicesC(ED), eigen_vals = eigen(bit_rates)$values, eigen_vecs = eigen(bit_rates)$vectors, inverse_vecs = solve(eigen_vecs)){
  st_root <- unname(st_labels[1,1]) #st_labels[1,] always corresponds to st_root by construction
  st_leaves <- st_labels[!(st_labels[,8] %in% st_labels[,1]),1]

  #Identify migrations with child and parent inside subtree
  rm_mig_pos <- (is.na(ED[,9])) & (ED[,8] %in% st_labels[,1]) & (ED[,7] %in% st_labels[,1])
  ED_rm_rows <- NI[ED[rm_mig_pos,1]]

  internal_coal_nodes <- st_labels[!(st_labels[,1] %in% c(st_root, st_leaves)),1]

  if (is.na(st_labels[1,7])){ #If st_root is global root and not fully determined
    if (all(st_labels[1, 8:9] %in% st_labels[,1])){
      internal_coal_nodes <- append(internal_coal_nodes, st_root)
    }
  }

  int_coal_node_row <- NI[internal_coal_nodes]

  #### Sample new deme at internal coalescent nodes
  n_deme <- nrow(bit_rates)
  ED[int_coal_node_row, 5] <- sample(1:n_deme, length(internal_coal_nodes), TRUE)

  ### Complete migration history
  max_label <- ED[nrow(ED), 1] #max(ED[,1]) is label of final row
  log_like <- 0 #Proposal log likelihood

  for (row_id in 2 : nrow(st_labels)){ #Row 1 = st_root
    node_row <- NI[st_labels[row_id, 1]]
    node_deme <- ED[node_row, 5]

    parent_row <- NI[st_labels[row_id, 7]]
    parent_deme <- ED[parent_row, 5]

    edge_length <- ED[node_row, 6] - ED[parent_row, 6]

    trans_mat <- Re(eigen_vecs %*% diag(exp(eigen_vals * edge_length))  %*% inverse_vecs)

    mig_path <- ECctmc::sample_path_unif3(node_deme, parent_deme,
                                          0, edge_length,
                                          Q = bit_rates, P = trans_mat)



    ###### ECctmc::sample_path_unif returns t(mig_path) if 0 migrations added
    ###### BUG IN SOURCE CODE!!!
    if (nrow(mig_path) == ncol(mig_path)){
      if (mig_path[1,2] == edge_length){
        mig_path <- t(mig_path)
      }
    }

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

  return(list(proposal=ED, prop_prob = log_like))
}


#' MultiTypeTree Probability Density
#'
#' Computes the relevant probability density to compute transition kernels for
#' the MultiTypeTree node retype move
#'
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param bit_mig_mat Matrix of initial backwards-in-time migration rates for the MCMC chain
#' @param NI Vector of row indices corresponding to which row of ED corresponds to each node label
#' @param st_labels Reduced extended data structure consisting of the rows of ED corresponding to coalescent nodes in the selected subtree only (local likelihood only)
#'
#' @return List consisting of the log.likelihood and likelihood of the structured genealogy
#'
#' @export

MTT_proposal_like <- function(ED, bit_mig_mat, NI = NodeIndicesC(ED)){
  parent_rows <- NI[ED[,2]]
  edge_lengths <- ED[,6] - ED[parent_rows, 6]

  bit_rowsums <- rowSums(bit_mig_mat)
  log_like <- 0


  parent_demes <- ED[parent_rows, 5]
  for (i in 1 : nrow(ED)){
    if (!is.na(parent_rows[i])){
      current_deme <- ED[i,5]
      log_like <- log_like - edge_lengths[i] * bit_rowsums[current_deme]

      if (parent_demes[i] != current_deme){
        log_like <- log_like + log(bit_mig_mat[current_deme, parent_demes[i]])
      }
    }
  }

  non_mig_nodes <- which((!is.na(ED[,4])) | (is.na(ED[,3])))
  bit_rates <- bit_mig_mat
  diag(bit_rates) <- 0
  diag(bit_rates) <- - rowSums(bit_mig_mat)

  for (row_id in non_mig_nodes){
    parent_row <- NI[ED[row_id, 7]]
    if (!is.na(parent_row)){
      edge_length <- ED[row_id, 6] - ED[parent_row, 6]
      trans_mat <- expm::expm(bit_rates * edge_length) #Use eigen decomposition!!

      log_like <- log_like - log(trans_mat[ED[row_id, 5], ED[parent_row, 5]])
    }
  }

  return(list(log.likelihood = log_like, likelihood = exp(log_like)))
}

#' @rdname MTT_proposal_like
#' @export

MTT_proposal_like_eigen <- function(ED, bit_rates, NI = NodeIndicesC(ED), eigen_vals = eigen(bit_rates)$values, eigen_vecs = eigen(bit_rates)$vectors, inverse_vecs = solve(eigen_vecs)){
  parent_rows <- NI[ED[,2]]
  edge_lengths <- ED[,6] - ED[parent_rows, 6]
  log_like <- 0


  parent_demes <- ED[parent_rows, 5]
  for (i in 1 : nrow(ED)){
    if (!is.na(parent_rows[i])){
      current_deme <- ED[i,5]
      log_like <- log_like + edge_lengths[i] * bit_rates[current_deme, current_deme] #bit_rates[i,i] = rowSums(bit_mig_mat)[i]

      if (parent_demes[i] != current_deme){
        log_like <- log_like + log(bit_rates[current_deme, parent_demes[i]])
      }
    }
  }

  non_mig_nodes <- which((!is.na(ED[,4])) | (is.na(ED[,3])))

  for (row_id in non_mig_nodes){
    parent_row <- NI[ED[row_id, 7]]
    if (!is.na(parent_row)){
      edge_length <- ED[row_id, 6] - ED[parent_row, 6]
      trans_mat <- eigen_vecs %*% diag(exp(edge_length * eigen_vals)) %*% inverse_vecs
      log_like <- log_like - log(trans_mat[ED[row_id, 5], ED[parent_row, 5]])
    }
  }

  return(list(log.likelihood = log_like, likelihood = exp(log_like)))
}

#' @rdname MTT_proposal_like
#' @export

MTT_local_like <- function(ED, st_labels, bit_mig_mat, NI = NodeIndicesC(ED)){
  #Modify st_labels into valid ED structure (no reference below leaves or above root)
  st_labels <- ED[NI[st_labels[,1]],]
  st_labels[1,c(2,7)] <- NA #Parents of st_root (in row 1) = NA

  is_leaf <- is.na(st_labels[,3]) | #Leaf of ED
    !((st_labels[,8] %in% st_labels[,1]) | (st_labels[,9] %in% st_labels[,1])) #Neither child coalescent node is in st_labels
  st_labels[is_leaf, c(3:4, 8:9)] <- NA #Children of st_leaves = NA

  #Extract migration events in subtree and add to st_labels
  is_st_mig <- (ED[,7] %in% st_labels[,1]) & #Parent coalescent node in subtree
    ((ED[,8] %in% st_labels[,1]) | (ED[,9] %in% st_labels[,1])) & #Child coalescent node in subtree
    is.na(ED[,4]) #Is migration

  st_labels <- rbind(st_labels, ED[is_st_mig,])

  return(MTT_proposal_like(st_labels, bit_mig_mat, NodeIndicesC(st_labels)))
}

#' @rdname MTT_proposal_like
#' @export

MTT_local_like_eigen <- function(ED, st_labels, bit_rates, NI = NodeIndicesC(ED), eigen_vals = NULL, eigen_vecs = NULL, inverse_vecs = NULL){
  if (is.null(eigen_vals) || is.null(eigen_vecs)){
    eigen_decomp <- eigen(bit_rates)
    eigen_vals <- eigen_decomp$values
    eigen_vecs <- eigen_decomp$vectors
    inverse_vecs <- solve(eigen_vecs)
  }

  #Modify st_labels into valid ED structure (no reference below leaves or above root)
  st_labels <- ED[NI[st_labels[,1]],]
  st_labels[1,c(2,7)] <- NA #Parents of st_root (in row 1) = NA

  is_leaf <- is.na(st_labels[,3]) | #Leaf of ED
    !((st_labels[,8] %in% st_labels[,1]) | (st_labels[,9] %in% st_labels[,1])) #Neither child coalescent node is in st_labels
  st_labels[is_leaf, c(3:4, 8:9)] <- NA #Children of st_leaves = NA

  #Extract migration events in subtree and add to st_labels
  is_st_mig <- (ED[,7] %in% st_labels[,1]) & #Parent coalescent node in subtree
    ((ED[,8] %in% st_labels[,1]) | (ED[,9] %in% st_labels[,1])) & #Child coalescent node in subtree
    is.na(ED[,4]) #Is migration

  st_labels <- rbind(st_labels, ED[is_st_mig,])

  return(MTT_proposal_like_eigen(st_labels, bit_rates, NodeIndicesC(st_labels), eigen_vals, eigen_vecs, inverse_vecs))
}

#' MultiTypeTree MCMC
#'
#' Runs a MCMC chain drawing joint samples of migration histories, backwards-in-time migration rates and coalescent rates.
#' MCMC chain runs with the node retype operator from MultiTypeTree for migration history updates and Gamma Gibbs moves for migration and coalescent rates
#'
#' @param N Number of MCMC iterations
#' @param ED Extended data representation of phylogenetic tree and initial migration history
#' @param coal_rate Vector of initial coalescent rates for the MCMC chain
#' @param bit_mig_mat Matrix of initial backwards-in-time migration rates for the MCMC chain
#' @param st_depth Number of generations below a coalescent node to include in the subtree for the node retype move
#' @param cr_mode Mode for prior on coalescent rates
#' @param cr_var Variance for prior on coalescent rates
#' @param mm_mode Mode for prior on backwards-in-time migration rates
#' @param mm_var Variance for prior on backwards-in-time migration rates
#' @param thin Thinning increment for the MCMC chain. A sample is saved every \texttt{thin} iterations
#' @param prop_rates Relative rates for proposals of migration history updates to coalescent rate updates to backwards-in-time migration rates updates
#' @param output_dir Directory to output files to for (thinned) results from the MCMC chain
#' @param run_name Name for run in output_dir
#'
#' @return Output files output_dir/run_name.freq, output_dir/run_name.log and output_dir/run_name.trees giving the acceptance rates of each MCMC move, thinned samples of coalescent and migration rates, and thinned samples of trees respectively
#'
#' @export

MTT_node_retype_MCMC <- function(N = 1e6,
                                 ED, coal_rate, bit_mig_mat,
                                 st_depth = 1, cr_mode = 1, cr_var = 1, mm_mode = 0.05, mm_var = 0.5,
                                 thin = 1e3, prop_rates = c(100, 1, 1),
                                 output_dir = '~', run_name = 'MTT_R'){
  pb <- txtProgressBar(0, N, style = 3)

  # Prior parameters
  cr_rate <- (cr_mode + sqrt(cr_mode^2 + 4 * cr_var))/(2 * cr_var)
  cr_shape <- 1 + cr_mode * cr_rate
  mm_rate <- (mm_mode + sqrt(mm_mode^2 + 4 * mm_var))/(2 * mm_var)
  mm_shape <- 1 + mm_mode * mm_rate

  bit_rates <- bit_mig_mat
  diag(bit_rates) <- - rowSums(bit_mig_mat)

  # Initial likelihoods and priors
  ED_NI <- NodeIndicesC(ED)
  ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)

  mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
  cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))

  # Initialise output files
  freq_file <- file.path(output_dir,
                         paste0(run_name, ".freq"))
  freq <- matrix(0, 2, 3,
                 dimnames = list(c("#accept", "#reject"), c("Node_Retype", "CR", "MM")))

  # log_file stores mig rates and coal rates etc
  log_file <- file.path(output_dir,
                        paste0(run_name, ".log"))
  # Column names
  cat("sample",
      "likelihood",
      "posterior",
      paste("coal_rate_", 1:n_deme, sep = "", collapse =","),
      paste("backward_migration_rate_", as.vector(outer(1:n_deme, 1:n_deme, paste, sep = "_")[-(1 + 0:(n_deme - 1) * (n_deme + 1))]), sep = ""),
      file = log_file, sep =",")

  # First row
  cat(paste0("\n", 0), #sample
      ED_SC, # likelihood
      ED_SC + mm_prior + cr_prior, # posterior
      coal_rate, #coal_rate
      as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
      file = log_file, sep =",", append = TRUE)

  # tree_file stores posterior sampled trees
  tree_file <- file.path(output_dir,
                         paste0(run_name, ".trees"))
  phylo <- ed.to.phylo(ED)
  treedata <- treeio::as.treedata(phylo)
  treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))

  # .nexus file header
  header <- capture.output(write.nexus(phylo, file = stdout(), translate = TRUE))
  header[2] <- paste("[R-package scoal, ", date(), "]\n\n", sep = "")
  cat(header[-(length(header)-0:1)], file = tree_file, sep = "\n")

  #Initial tree
  cat("\tTREE STATE_0 = ",
      treeio::write.beast.newick(treedata), "\n",
      file = tree_file, append = TRUE, sep = "")

  for (x in 1 : N){
    move_id <- sample(1:3, 1, prob = prop_rates)

    if (move_id == 1){
      subtree <- MTT_st_coal_node(ED, st_depth)
      proposal <- MTT_node_retype(ED, subtree$st_labels, bit_rates, ED_NI)

      # Early acceptance if prop == ED
      if ((nrow(ED) == nrow(proposal$proposal)) && (all(na.omit(as.vector(ED == proposal$proposal))))){
        freq[1, move_id] <- freq[1, move_id] + 1
      } else {
        prop_NI <- NodeIndicesC(proposal$proposal)
        prop_SC <- SC_like_C(proposal$proposal, coal_rate, bit_mig_mat, prop_NI)
        log_AR <- min(0, prop_SC - ED_SC +
                      MTT_local_like(ED, subtree$st_labels, bit_mig_mat, ED_NI)$log.likelihood -
                      # proposal$prop_prob
                      MTT_local_like(proposal$proposal, subtree$st_labels, bit_mig_mat, prop_NI)$log.likelihood
        )

        if (log(runif(1)) < log_AR){
          # Accept
          freq[1, move_id] <- freq[1, move_id] + 1

          ED <- proposal$proposal
          ED_NI <- prop_NI
          ED_SC <- prop_SC
        }
      }
    } else {
      if (move_id == 2){
        coal_rate <- scaled_coal_rate_gibbs_update(ED, 1, n_deme, ED_NI, cr_shape, cr_rate)
        cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))
      } else if (move_id == 3){
        bit_mig_mat <- scaled_mig_mat_gibbs_update(ED, 1, n_deme, ED_NI, mm_shape, mm_rate)

        bit_rates <- bit_mig_mat
        diag(bit_rates) <- - rowSums(bit_mig_mat)
        mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
      }

      ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)
      freq[1, move_id] <- freq[1, move_id] + 1
    }

    freq[2, move_id] <- freq[2, move_id] + 1

    setTxtProgressBar(pb, x)
    if (x %% thin == 0){ #x (mod thin) = 0, i.e. thin-many iterations have passed
      cat(paste0("\n", x), #sample
          ED_SC, # likelihood
          ED_SC + mm_prior + cr_prior, # posterior
          coal_rate, #coal_rate
          as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
          file = log_file, sep =",", append = TRUE)

      phylo <- ed.to.phylo(ED)
      treedata <- treeio::as.treedata(phylo)
      treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))
      cat("\tTREE STATE_", x, " = ",
          treeio::write.beast.newick(treedata), "\n",
          file = tree_file, append = TRUE, sep = "")

      write.table(freq,
                  file = freq_file,
                  row.names = c('#ACCEPT', '#TOTAL'),
                  col.names = TRUE)
    }

  }
  close(pb)
  cat("END;",
      file = tree_file, append = TRUE, sep = "")
}

#' @rdname MTT_node_retype_MCMC
#' @export

MTT_node_retype_MCMC_eigen <- function(N = 1e6,
                                       ED, coal_rate, bit_mig_mat,
                                       st_depth = 1, cr_mode = 1, cr_var = 1, mm_mode = 0.05, mm_var = 0.5,
                                       thin = 1e3, tree_thin=1e3,
                                       prop_rates = c(100, 1, 1),
                                       output_dir = '~', run_name = 'MTT_R'){
  # Prior parameters
  cr_rate <- (cr_mode + sqrt(cr_mode^2 + 4 * cr_var))/(2 * cr_var)
  cr_shape <- 1 + cr_mode * cr_rate
  mm_rate <- (mm_mode + sqrt(mm_mode^2 + 4 * mm_var))/(2 * mm_var)
  mm_shape <- 1 + mm_mode * mm_rate

  bit_rates <- bit_mig_mat
  diag(bit_rates) <- - rowSums(bit_mig_mat)
  n_deme <- nrow(bit_rates)

  eigen_decomp <- eigen(bit_rates)
  eigen_vecs <- eigen_decomp$vectors
  eigen_vals <- eigen_decomp$values
  inverse_vecs <- solve(eigen_vecs)

  # Initial likelihoods and priors
  ED_NI <- NodeIndicesC(ED)
  ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)

  mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
  cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))

  # Initialise output files
  freq_file <- file.path(output_dir,
                         paste0(run_name, ".freq"))
  freq <- matrix(0, 2, 3,
                 dimnames = list(c("#accept", "#reject"), c("Node_Retype", "CR", "MM")))

  # log_file stores mig rates and coal rates etc
  log_file <- file.path(output_dir,
                        paste0(run_name, ".log"))
  # Column names
  cat("sample",
      "likelihood",
      "posterior",
      paste("coal_rate_", 1:n_deme, sep = "", collapse =","),
      paste("backward_migration_rate_", as.vector(outer(1:n_deme, 1:n_deme, paste, sep = "_")[-(1 + 0:(n_deme - 1) * (n_deme + 1))]), sep = ""),
      file = log_file, sep =",")

  # First row
  cat(paste0("\n", 0), #sample
      ED_SC, # likelihood
      ED_SC + mm_prior + cr_prior, # posterior
      coal_rate, #coal_rate
      as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
      file = log_file, sep =",", append = TRUE)

  # tree_file stores posterior sampled trees
  tree_file <- file.path(output_dir,
                         paste0(run_name, ".trees"))
  phylo <- ed.to.phylo(ED)
  treedata <- treeio::as.treedata(phylo)
  treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))

  # .nexus file header
  header <- capture.output(write.nexus(phylo, file = stdout(), translate = TRUE))
  header[2] <- paste("[R-package scoal, ", date(), "]\n\n", sep = "")
  cat(header[-(length(header)-0:1)], file = tree_file, sep = "\n")

  #Initial tree
  cat("\tTREE STATE_0 = ",
      treeio::write.beast.newick(treedata), "\n",
      file = tree_file, append = TRUE, sep = "")

  for (x in 1 : N){
    move_id <- sample(1:3, 1, prob = prop_rates)

    if (move_id == 1){
      subtree <- MTT_st_coal_node(ED, st_depth)
      proposal <- MTT_node_retype_eigen(ED, subtree$st_labels, bit_rates, ED_NI, eigen_vals, eigen_vecs, inverse_vecs)

      # Early acceptance if prop == ED
      if ((nrow(ED) == nrow(proposal$proposal)) && (all(na.omit(as.vector(ED == proposal$proposal))))){
        freq[1, move_id] <- freq[1, move_id] + 1
      } else {
        prop_NI <- NodeIndicesC(proposal$proposal)
        prop_SC <- SC_like_C(proposal$proposal, coal_rate, bit_mig_mat, prop_NI)
        log_AR <- min(0, Re(prop_SC - ED_SC +
                              MTT_local_like_eigen(ED, subtree$st_labels, bit_rates, ED_NI, eigen_vals, eigen_vecs, inverse_vecs)$log.likelihood -
                              # proposal$prop_prob
                              MTT_local_like_eigen(proposal$proposal, subtree$st_labels, bit_rates, prop_NI, eigen_vals, eigen_vecs, inverse_vecs)$log.likelihood
        ))

        if (log(runif(1)) < log_AR){
          # Accept
          freq[1, move_id] <- freq[1, move_id] + 1

          ED <- proposal$proposal
          ED_NI <- prop_NI
          ED_SC <- prop_SC
        }
      }
    } else {
      if (move_id == 2){
        coal_rate <- scaled_coal_rate_gibbs_update(ED, 1, n_deme, ED_NI, cr_shape, cr_rate)
        cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))
      } else if (move_id == 3){
        bit_mig_mat <- scaled_mig_mat_gibbs_update(ED, 1, n_deme, ED_NI, mm_shape, mm_rate)

        bit_rates <- bit_mig_mat
        diag(bit_rates) <- - rowSums(bit_mig_mat)

        eigen_decomp <- eigen(bit_rates)
        eigen_vecs <- eigen_decomp$vectors
        eigen_vals <- eigen_decomp$values
        inverse_vecs <- solve(eigen_vecs)

        mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
      }

      ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)
      freq[1, move_id] <- freq[1, move_id] + 1
    }

    freq[2, move_id] <- freq[2, move_id] + 1


    if (x %% thin == 0){
      # Write iteration, likelihood and posterior to stdout()
      cat(x, # Sample
          sprintf('%.03f', ED_SC), # Likelihood
          sprintf('%.03f', ED_SC + mm_prior + cr_prior), # posterior
          sep = '\t')

      # Write continuous parameters to .log file
      cat(paste0("\n", x), #sample
          ED_SC, # likelihood
          ED_SC + mm_prior + cr_prior, # posterior
          coal_rate, #coal_rate
          as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
          file = log_file, sep =",", append = TRUE)

      # Update .freq file (Overwrites existing file entirely)
      write.table(freq,
                  file = freq_file,
                  row.names = c('#ACCEPT', '#TOTAL'),
                  col.names = TRUE)
    }

    if (x %% tree_thin == 0){
      # Write current tree to .trees file
      phylo <- ed.to.phylo(ED)
      treedata <- treeio::as.treedata(phylo)
      treedata@data <- tidytree::tibble(type = paste0("\"", phylo$node.deme, "\""), node = 1:length(phylo$node.deme))
      cat("\tTREE STATE_", x, " = ",
          treeio::write.beast.newick(treedata), "\n",
          file = tree_file, append = TRUE, sep = "")
    }
  }
  cat("END;",
      file = tree_file, append = TRUE, sep = "")
}
