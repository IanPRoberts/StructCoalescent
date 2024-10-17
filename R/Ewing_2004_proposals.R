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
#'
#' @export

migration_birth <- function(ED, n_deme, fix_leaf_deme = TRUE, node_indices){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  parent.rows <- node_indices[ED[,2]]
  parent.times <- ED[parent.rows, 6]
  parent.times[root.node] <- 0
  edge.length <- ED[,6] - parent.times

  tree.length <- sum(edge.length)  #Total tree length
  new.location <- runif(1, 0, tree.length)

  child.row <- min(which(cumsum(edge.length) >= new.location))
  child.node <- ED[child.row, 1]

  parent.node <- ED[child.row, 2]
  parent.row <- node_indices[parent.node]

  #Calculating subtree containing edge <new.node, child.node> terminating in migration or leaf nodes
  subtree.nodes <- child.node
  active.nodes <- child.node
  while (any(active.nodes %in% coalescence.nodes)){
    active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]

    active.rows <- node_indices[active.nodes]
    subtree.nodes <- c(subtree.nodes, ED[active.rows, 3:4]) #Add children of active.nodes to subtree
    active.nodes <- ED[active.rows, 3:4]
  }

  if ((fix_leaf_deme == TRUE) && (any(subtree.nodes %in% leaf.nodes))){
    # REJECT
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
  } else{
    #Continue proposal
    subtree.leaves <- subtree.nodes[subtree.nodes %in% migration.nodes]
    old.deme <- ED[child.row, 5]

    proposal.deme <- sample.vector((1:n_deme)[- old.deme], 1)  #Propose deme update without accounting for surrounding demes

    subtree.leaf.rows <- node_indices[subtree.leaves]
    subtree.leaf.children <- ED[subtree.leaf.rows, 3]
    subtree.leaf.children.rows <- node_indices[subtree.leaf.children]

    if (any(ED[subtree.leaf.children.rows, 5] == proposal.deme)){
      #REJECT
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
    }

    #Update deme across subtree
    subtree.rows <- node_indices[subtree.nodes]
    ED[subtree.rows, 5] <- proposal.deme

    new.node <- max(ED[,1]) + 1

    #Update child.row with new parent
    ED[child.row, 2] <- new.node

    #Update parent row with new child
    which.child <- which(ED[parent.row, 3:4] == child.node)
    ED[parent.row, 2 + which.child] <- new.node  #which.child either 1 or 2 so add 2 to get to col no.

    #Add row for new.node
    cumulative.edge.length <- c(0, cumsum(edge.length)[-length(edge.length)])  #Offset cumsum(edge.length) by 1 in case new.location is on first edge in edge length
    new.age <- ED[child.row, 6] - (new.location - max(cumulative.edge.length[cumulative.edge.length < new.location])) #Age of new.node (currently have distance along tree from new.location)

    ED <- rbind(ED, c(new.node, parent.node, child.node, NA, old.deme, new.age))
    prop.ratio <- (n_deme - 1) * tree.length / (length(migration.nodes) + 1)
    log.prop.ratio <- log(n_deme - 1) + log(tree.length) - log(length(migration.nodes) + 1)

    if (new.node > length(node_indices)){
      new.node_indices <- numeric(new.node)
      new.node_indices[1:length(node_indices)] <- node_indices
    } else{
      new.node_indices <- node_indices
    }
    new.node_indices[new.node] <- dim(ED)[1]

    return(list(ED= ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node_indices = new.node_indices))
  }
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
#'
#' @export

ed.mig.death <- function(ED, n_deme, fix_leaf_deme = TRUE, node_indices){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]
  coalescence.nodes <- coalescence.nodes[coalescence.nodes != root.node]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  parent.rows <- node_indices[ED[,2]]
  parent.times <- ED[parent.rows, 6]
  parent.times[root.node] <- 0
  edge.length <- ED[,6] - parent.times

  tree.length <- sum(edge.length)  #Total tree length

  if (length(migration.nodes) == 0){
    #REJECT
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
  }
  selected.node <- sample.vector(migration.nodes, 1)
  selected.row <- node_indices[selected.node]

  child.node <- ED[selected.row, 3]
  child.row <- node_indices[child.node]

  #Calculating subtree containing edge <new.node, child.node> terminating in migration or leaf nodes
  subtree.nodes <- child.node
  active.nodes <- child.node
  while (any(active.nodes %in% coalescence.nodes)){
    active.nodes <- active.nodes[active.nodes %in% coalescence.nodes]
    active.rows <- node_indices[active.nodes]
    subtree.nodes <- c(subtree.nodes, ED[active.rows, 3:4]) #Add children of active.nodes to subtree
    active.nodes <- ED[active.rows, 3:4]
  }

  if ((fix_leaf_deme == TRUE) && (any(subtree.nodes %in% leaf.nodes))){
    # REJECT
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
  } else{
    #Continue proposal
    subtree.leaves <- subtree.nodes[subtree.nodes %in% migration.nodes]

    subtree.leaf.rows <- node_indices[subtree.leaves]
    subtree.leaf.children <- ED[subtree.leaf.rows, 3]
    subtree.leaf.children.rows <- node_indices[subtree.leaf.children]
    forbidden.demes <- unique(ED[subtree.leaf.children.rows, 5])

    proposal.deme <- ED[selected.row, 5]  #Propose deme immediately above subtree
    if (proposal.deme %in% forbidden.demes){
      #REJECT - cannot change deme
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
    }

    parent.node <- ED[selected.row, 2]
    parent.row <- node_indices[parent.node]

    #Update deme across subtree
    subtree.rows <- node_indices[subtree.nodes]
    ED[subtree.rows, 5] <- proposal.deme

    ED[child.row, 2] <- parent.node  #Update child.node parent

    #Update parent.node child
    which.child <- which(ED[parent.row, 3:4] == selected.node)
    ED[parent.row, 2 + which.child] <- child.node

    node_indices[selected.node] <- 0
    index.changes <- (node_indices > selected.row)
    node_indices[index.changes] <- node_indices[index.changes] - 1

    ED <- ED[-selected.row,]
    prop.ratio <- length(migration.nodes)/ ((n_deme - 1) * tree.length)
    log.prop.ratio <- log(length(migration.nodes)) - log(n_deme - 1) - log(tree.length)
    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node_indices = node_indices))
  }
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
#'
#' @export

ed.mig.pair.birth <- function(ED, n_deme, node_indices){
  root.node <- which(is.na(ED[,2]))

  #Sample non-root node to obtain edge <sampled.node, node.parent>
  selected.node <- sample(ED[-root.node, 1], 1)
  selected.row <- node_indices[selected.node]
  parent.node <- ED[selected.row, 2]
  parent.row <- node_indices[parent.node]

  new.nodes <- max(ED[,1]) + c(1,2)  #New node IDs
  new.node.ages <- sort(runif(2, min = ED[parent.row, 6], max = ED[selected.row, 6]), decreasing = TRUE) #New node ages

  new.deme <- sample.vector((1:n_deme)[-ED[selected.row,5]], 1)

  ED[selected.row, 2] <- new.nodes[1]
  ED[parent.row, 2 + which(ED[parent.row,3:4] == selected.node)] <- new.nodes[2]
  ED <- rbind(ED,
              c(new.nodes[1], new.nodes[2], selected.node, NA, new.deme, new.node.ages[1]),
              c(new.nodes[2], parent.node, new.nodes[1], NA, ED[selected.row, 5], new.node.ages[2])
  )

  n.nodes <- dim(ED)[1]
  prop.ratio <- (n_deme - 1) * (n.nodes - 1) * (ED[selected.row, 6] - ED[parent.row, 6])^2 / (2 * (n.nodes + 1))
  log.prop.ratio <- log(n_deme - 1) + log(n.nodes - 1) + 2 * log(abs(ED[selected.row, 6] - ED[parent.row, 6])) - log(2) - log(n.nodes + 1)

  if (new.nodes[2] > length(node_indices)){
    new.node_indices <- numeric(new.nodes[2])
    new.node_indices[1:length(node_indices)] <- node_indices
  } else{
    new.node_indices <- node_indices
  }
  new.node_indices[new.nodes] <- c(dim(ED)[1] - 1, dim(ED)[1])
  return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node_indices = new.node_indices))
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
#'
#' @export

ed.mig.pair.death <- function(ED, n_deme, node_indices){
  root.node <- ED[is.na(ED[,2]), 1]

  #Sample non-root node to obtain edge <sampled.node, node.parent>
  selected.node <- sample(ED[-root.node, 1], 1)
  selected.row <- node_indices[selected.node]
  parent.node <- ED[selected.row, 2]
  parent.row <- node_indices[parent.node]

  if (((is.na(ED[selected.row, 4]) & (!is.na(ED[selected.row, 3])))) &&
      ((is.na(ED[parent.row, 4]) & (!is.na(ED[parent.row, 3])))) &&
      (ED[node_indices[ED[selected.row, 3]], 5] == ED[parent.row, 5])){
    #Both ends of the edge are migration nodes, and the demes are consistent to remove the pair of nodes
    parent2.node <- ED[parent.row, 2]
    parent2.row <- node_indices[parent2.node]
    child.node <- ED[selected.row, 3]
    child.row <- node_indices[child.node]

    n.nodes <- dim(ED)[1]
    prop.ratio <- 2 * (n.nodes - 1) / ((n.nodes - 3) * (n_deme - 1) * (ED[child.row, 6] - ED[parent2.row, 6])^2)
    log.prop.ratio <- log(2) + log(n.nodes - 1) - log(n.nodes - 3) - log(n_deme - 1) - 2 * log(abs(ED[child.row, 6] - ED[parent2.row, 6]))

    ED[parent2.row, 2 + which(ED[parent2.row, 3:4] == parent.node)] <- child.node
    ED[child.row, 2] <- parent2.node
    ED <- ED[-c(selected.row, parent.row),]

    node_indices[c(selected.node, parent.node)] <- 0
    index.changes.1 <- (node_indices > selected.row)
    index.changes.2 <- (node_indices > parent.row)
    node_indices[index.changes.1] <- node_indices[index.changes.1] - 1
    node_indices[index.changes.2] <- node_indices[index.changes.2] - 1

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node_indices = node_indices))
  } else {
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
  }
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
#'
#' @export

ed.coal.merge <- function(ED, n_deme, node_indices){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]

  selected.node <- sample.vector(coalescence.nodes, 1)
  selected.row <- node_indices[selected.node]

  child.nodes <- ED[selected.row, 3:4]
  child.rows <- node_indices[child.nodes]
  if (!all((is.na(ED[child.rows,4])) & (!is.na(ED[child.rows,3])))){
    #REJECT as cannot delete coalescent node
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
  }

  child.child <- ED[child.rows, 3]
  child.child.rows <- node_indices[child.child]
  exterior.demes <- ED[child.child.rows, 5]

  if (exterior.demes[1] != exterior.demes[2]){
    #REJECT as exterior demes are not the same
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
  }

  if (selected.node == root.node){
    prop.ratio <- 1 / ((n_deme - 1)* (ED[child.child.rows[1], 6] - ED[selected.row, 6]) * (ED[child.child.rows[2], 6] - ED[selected.row, 6]) )
    log.prop.ratio <- -log(n_deme - 1) - log(abs(ED[child.child.rows[1], 6] - ED[selected.row, 6])) - log(abs(ED[child.child.rows[2], 6] - ED[selected.row, 6]))

    ED[selected.row, c(3:5)] <- c(child.child, exterior.demes[1])
    ED[child.child.rows, 2] <- selected.node
    ED <- ED[!ED[,1] %in% child.nodes, ] #Remove child migration nodes

    node_indices[child.nodes] <- 0
    index.changes.1 <- (node_indices > child.rows[1])
    index.changes.2 <- (node_indices > child.rows[2])
    node_indices[index.changes.1] <- node_indices[index.changes.1] - 1
    node_indices[index.changes.2] <- node_indices[index.changes.2] - 1

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node_indices = node_indices))
  } else{
    parent.node <- ED[selected.row, 2]
    parent.row <- node_indices[parent.node]

    new.node <- max(ED[,1]) + 1
    new.node.time <- runif(1, ED[parent.row, 6], ED[selected.row, 6])
    old.deme <- ED[selected.row, 5]

    ED[selected.row, c(2:5)] <- c(new.node, child.child, exterior.demes[1])
    ED[c(child.child.rows[1], child.child.rows[2]), 2] <- selected.node

    which.child <- which(ED[parent.row, 3:4] == selected.node)
    ED[parent.row, 2 + which.child] <- new.node

    prop.ratio <- (ED[selected.row, 6] - ED[parent.row ,6]) / ((ED[child.child.rows[1], 6] - ED[selected.row, 6]) * (ED[child.child.rows[2], 6] - ED[selected.row, 6]))
    log.prop.ratio <- log(abs(ED[selected.row, 6] - ED[parent.row ,6])) - log(abs(ED[child.child.rows[1], 6] - ED[selected.row, 6])) - log(abs(ED[child.child.rows[2], 6] - ED[selected.row, 6]))


    ED <- ED[!ED[,1] %in% child.nodes, ] #Remove child migration nodes
    ED <- rbind(ED, c(new.node, parent.node, selected.node, NA, old.deme, new.node.time)) #Add new parent migration node

    node_indices[child.nodes] <- 0
    index.changes.1 <- (node_indices > child.rows[1])
    index.changes.2 <- (node_indices > child.rows[2])
    node_indices[index.changes.1] <- node_indices[index.changes.1] - 1
    node_indices[index.changes.2] <- node_indices[index.changes.2] - 1

    if (new.node > length(node_indices)){
      new.node_indices <- numeric(new.node)
      new.node_indices[1:length(node_indices)] <- node_indices
    } else{
      new.node_indices <- node_indices
    }
    new.node_indices[new.node] <- dim(ED)[1]

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node_indices = new.node_indices))
  }
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
#'
#' @export

ed.coal.split <- function(ED, n_deme, node_indices){
  root.node <- ED[is.na(ED[,2]), 1]
  coalescence.nodes <- ED[!is.na(ED[,4]),1]

  selected.node <- sample.vector(coalescence.nodes, 1)
  selected.row <- node_indices[selected.node]

  if (selected.node == root.node){
    child.nodes <- ED[selected.row, 3:4]
    child.rows <- node_indices[child.nodes]

    new.nodes <- max(ED[,1]) + c(1,2)
    new.node.times <- sapply(ED[child.rows, 6], runif, n = 1, min = ED[selected.row, 6])

    proposal.deme <- sample.vector((1:n_deme)[- ED[selected.row, 5]], 1)  #Sample new deme not equal to current deme of the root

    prop.ratio <- (n_deme - 1)* (ED[child.rows[1], 6] - ED[selected.row, 6]) * (ED[child.rows[2], 6] - ED[selected.row, 6])
    log.prop.ratio <- log(n_deme - 1) + log(abs(ED[child.rows[1], 6] - ED[selected.row, 6])) + log(abs(ED[child.rows[2], 6] - ED[selected.row, 6]))
    ED[selected.row, 3:5] <- c(new.nodes, proposal.deme)
    ED[child.rows, 2] <- new.nodes
    ED <- rbind(ED, matrix(c(new.nodes, rep(selected.node,2), child.nodes, rep(NA,2), rep(proposal.deme, 2), new.node.times), 2))

    if (new.nodes[2] > length(node_indices)){
      new.node_indices <- numeric(new.nodes[2])
      new.node_indices[1:length(node_indices)] <- node_indices
    } else{
      new.node_indices <- node_indices
    }
    new.node_indices[new.nodes] <- c(dim(ED)[1] - 1, dim(ED)[1])

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node_indices = new.node_indices))
  } else{
    parent.node <- ED[selected.row, 2]
    parent.row <- node_indices[parent.node]

    if (!((is.na(ED[parent.row,4])) & (!is.na(ED[parent.row,3])))){
      #REJECT
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
    }

    child.nodes <- ED[selected.row, 3:4]
    child.rows <- node_indices[child.nodes]

    parent.parent.node <- ED[parent.row, 2]
    parent.parent.row <- node_indices[parent.parent.node]
    new.nodes <- max(ED[,1]) + c(1,2)
    new.node.times <- sapply(ED[child.rows, 6], runif, n = 1, min = ED[selected.row, 6])

    prop.ratio <- ((ED[child.rows[1], 6] - ED[selected.row, 6]) * (ED[child.rows[2], 6] - ED[selected.row, 6])) / (ED[selected.row, 6] - ED[parent.parent.row ,6])
    log.prop.ratio <- log(abs(ED[child.rows[1], 6] - ED[selected.row, 6])) + log(abs(ED[child.rows[2], 6] - ED[selected.row, 6])) - log(abs(ED[selected.row, 6] - ED[parent.parent.row ,6]))

    which.child <- which(ED[parent.parent.row, 3:4] == parent.node)
    ED[parent.parent.row, 2 + which.child] <- selected.node
    ED[child.rows, 2] <- new.nodes
    ED[selected.row, 2:5] <- c(parent.parent.node, new.nodes, ED[parent.row, 5])
    ED <- ED[-parent.row, ]
    ED <- rbind(ED, matrix(c(new.nodes, rep(selected.node, 2), child.nodes, rep(NA, 2), rep(ED[selected.row, 5],2), new.node.times), 2))

    node_indices[parent.node] <- 0
    index.changes <- (node_indices > parent.row)
    node_indices[index.changes] <- node_indices[index.changes] - 1

    if (new.nodes[2] > length(node_indices)){
      new.node_indices <- numeric(new.nodes[2])
      new.node_indices[1:length(node_indices)] <- node_indices
    } else{
      new.node_indices <- node_indices
    }
    new.node_indices[new.nodes] <- c(dim(ED)[1] - 1, dim(ED)[1])

    return(list(ED = ED, prop.ratio = prop.ratio, log.prop.ratio = log.prop.ratio, node_indices = new.node_indices))
  }
}

#' Block Recolouring MCMC Proposal
#'
#' Performs a block recolouring MCMC proposal on a structured coalescent
#' migration history. A block of the migration history currently in a single
#' deme is migrated into a new deme selected uniformly from all other demes. If
#' an inconsistent migration history is obtained, the move will be rejected.
#'
#' @param ED Extended data object; matrix with columns Node ID, parent, child 1, child 2, deme, node age
#' @param n_deme Number of distinct demes in the population
#' @param fix_leaf_deme Logical; if TRUE, the deme of each leaf in the tree is fixed. Any proposal which attempts to change a leaf deme is automatically rejected
#'
#' @return Updated extended data object with the proposal from the block recolouring proposal
#'
#' @export

ed.block.recolour <- function(ED, n_deme, fix_leaf_deme = TRUE, node_indices){
  root.node <- ED[is.na(ED[,2]), 1]
  leaf.nodes <- ED[(is.na(ED[,3])) & (is.na(ED[,4])), 1]
  migration.nodes <- ED[ is.na(ED[,4]) & (!is.na(ED[,3])) ,1]

  if (length(migration.nodes) == 0){
    if (fix_leaf_deme == TRUE){
      #REJECT - no migration nodes and fixed leaf demes -> cannot change deme
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
    } else{ #No migration nodes -> all same deme
      proposal.deme <- sample.vector((1:n_deme)[-ED[1, 5]], 1)
      ED[,5] <- proposal.deme
      return(list(ED = ED, prop.ratio = 1, log.prop.ratio = 0, node_indices = node_indices))
    }
  }

  parent.rows <- node_indices[ED[,2]]
  parent.times <- ED[parent.rows, 6]
  parent.times[root.node] <- 0
  edge.length <- ED[,6] - parent.times

  tree.length <- sum(edge.length)  #Total tree length
  new.location <- runif(1, 0, tree.length)  #Uniform location along tree.length

  child.row <- min(which(cumsum(edge.length) >= new.location))
  child.node <- ED[child.row, 1]

  #Sweep upwards until hitting the root or a migration node
  subtree.root <- ED[child.row, 2]
  while (! any(subtree.root == c(migration.nodes, root.node))){
    subtree.root.row <- node_indices[subtree.root]
    subtree.root <- ED[subtree.root.row, 2]
  }

  #Sweep downwards until hitting migration nodes
  subtree.root.row <- node_indices[subtree.root]
  if (subtree.root == root.node){
    active.nodes <- ED[subtree.root.row, 3:4]
    subtree.nodes <- c(subtree.root, active.nodes)
  } else{
    active.nodes <- ED[subtree.root.row, 3]
    subtree.nodes <- c(subtree.root, active.nodes)
  }

  while (! all(active.nodes %in% migration.nodes)){
    active.nodes <- active.nodes[! (active.nodes %in% c(migration.nodes, leaf.nodes))]
    active.rows <- node_indices[active.nodes]
    subtree.nodes <- c(subtree.nodes, ED[active.rows, 3:4])
    active.nodes <- ED[active.rows, 3:4]
  }

  if ((fix_leaf_deme == TRUE) && (any(subtree.nodes %in% leaf.nodes))){
    # REJECT
    return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
  } else{
    #Continue proposal
    subtree.leaves <- subtree.nodes[subtree.nodes %in% migration.nodes]
    if (subtree.root %in% subtree.leaves){
      subtree.leaves <- subtree.leaves[subtree.leaves != subtree.root]
    }
    old.deme <- ED[child.row, 5]

    proposal.deme <- sample.vector((1:n_deme)[- old.deme], 1)  #Propose deme update without accounting for surrounding demes

    subtree.leaf.rows <- node_indices[subtree.leaves]
    subtree.leaf.children <- ED[subtree.leaf.rows, 3]
    subtree.leaf.children.rows <- node_indices[subtree.leaf.children]

    if (any(ED[subtree.leaf.children.rows, 5] == proposal.deme)){
      #REJECT - induces self-migration
      return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
    }

    if (subtree.root != root.node){
      if (ED[subtree.root.row, 5] == proposal.deme){
        #REJECT
        return(list(ED = ED, prop.ratio = 0, log.prop.ratio = -Inf, node_indices = node_indices))
      }
    }

    #Update deme across subtree
    if (subtree.root == root.node){
      subtree.rows <- node_indices[subtree.nodes]
    } else{
      subtree.rows <- node_indices[subtree.nodes[subtree.nodes != subtree.root]]
    }
    ED[subtree.rows, 5] <- proposal.deme

    return(list(ED = ED, prop.ratio = 1, log.prop.ratio = 0, node_indices = node_indices))
  }
}
