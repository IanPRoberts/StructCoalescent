#' Simulation of Heterochronous Structured Phylogenies
#'
#' Simulates a dated phylogenetic tree under the structured coalescent model
#'
#' @param n Number of samples
#' @param d Number of demes
#' @param tip_data nx3 matrix with first column giving tip labels, second column giving tip ages and third column giving tip demes
#' @param coalescent_rates Numeric vector giving coalescent rates (inverse effective population sizes) of each deme
#' @param migration_matrix Numeric matrix giving backwards-in-t migration rates between each pair of demes
#'
#' @return An object of class '\code{strphylo}'
#'
#' @export

rstrphylo <- function(n, d, tip_data=cbind(1:n,0,rep_len(1:d,n)), coalescent_rates=1/effective_population_sizes, migration_matrix, effective_population_sizes=1/coalescent_rates){
  #Order by decreasing tip age
  tip_order <- order(tip_data[,2], decreasing = TRUE)
  tip_data <- tip_data[tip_order,]

  ED <- cbind('Node_ID'=1:n,
              'Parent'=NA,
              'Child_1'=NA,
              'Child_2'=NA,
              'Deme'=tip_data[,3],
              'Node_age'=tip_data[,2],
              'Parent_coal'=NA,
              'Child_coal_1'=NA,
              'Child_coal_2'=NA)

  t <- ED[1, 6]
  node_ID <- n+1 #Next node ID, reserving n+1 for root

  active_rows <- which(ED[,6] >= t)
  next_tip_age <- max(ED[ED[,6] < t, 6], -Inf)
  count <- 1

  k <- sapply(1:d, function(x) sum(ED[active_rows,5] == x)) #Number of active lineages in each deme

  while ((next_tip_age > -Inf) | (length(active_rows) > 1)){
    total_migration_rate <- sum(migration_matrix * k)
    total_coalescent_rate <- sum(k * (k-1) * coalescent_rates) / 2
    event_rate <- total_migration_rate + total_coalescent_rate

    sampling_prob <- 1 - pexp(t - next_tip_age, rate=event_rate) #Probability next event adds more tips

    if (runif(1) < sampling_prob){ #Sampling event
      t <- next_tip_age
      active_rows <- append(active_rows,
                            which(ED[1:n,6] == t)) #Add leaves at time t to active_rows
      next_tip_age <- max(ED[ED[,6] < t, 6], -Inf)
      k <- sapply(1:d, function(x) sum(ED[active_rows,5] == x)) #Number of active lineages in each deme
    }
    else{  #Event occurs before new nodes added
      t <- t - rexp_trunc(1, t - next_tip_age, event_rate)

      if (runif(1) <= total_coalescent_rate / event_rate){  #Coalescence event
        coalescence_deme <- sample.int(d, 1, prob = k * (k-1) * coalescent_rates / 2)
        ED <- rbind(ED,
                    c(node_ID,
                      NA,
                      sample(ED[active_rows[ED[active_rows, 5] == coalescence_deme], 1], 2),
                      coalescence_deme,
                      t,
                      rep(NA, 3)))
        k[coalescence_deme] <- k[coalescence_deme] - 1  #Number of lineages in coalescence deme decreases by 1
      }
      else{  #Migration event
        origin <- sample.int(d, 1, prob=k * rowSums(migration_matrix))
        target <- sample.int(d, 1, prob=migration_matrix[origin,])

        if (k[origin] == 1){
          ED <- rbind(ED,
                      c(node_ID,
                        NA,
                        ED[active_rows[ED[active_rows, 5] == origin],1],
                        NA,
                        target,
                        t,
                        rep(NA, 3)))
        } else {
          ED <- rbind(ED,
                      c(node_ID,
                        NA,
                        sample(ED[active_rows[ED[active_rows, 5] == origin],1], 1),
                        NA,
                        target,
                        t,
                        rep(NA, 3)))
        }

        k[origin] <- k[origin] - 1
        k[target] <- k[target] + 1
      }
      active_rows <- c(active_rows[!(ED[active_rows, 1] %in% ED[nrow(ED), 3:4])], node_ID)
      node_ID <- node_ID + 1
    }
  }

  #Add parents to ED
  for (row_id in 1 : nrow(ED)){
    for (child_id in 1:2){
      if (!is.na(ED[row_id, 2+child_id])){
        ED[ED[row_id, 2 + child_id], 2] <- ED[row_id, 1]
      }
    }
  }

  #Reorder ED into leaf -> root -> coal -> migrations
  ED <- ED[c(1:n,
             nrow(ED),
             which(!is.na(ED[-nrow(ED),4])),
             which((is.na(ED[,4])) & !is.na(ED[,3]))),]

  if (ED[n+1, 1] != n+1){
    max_label <- max(ED[,1])
    ED[,1:4][ED[,1:4] == n+1] <- max_label + 1
    ED[,1:4][ED[,1:4] == ED[n+1, 1]] <- n+1
  }

  class(ED) <- 'ED'
  strphylo <- as.strphylo(ED)
  strphylo$tip.label <- tip_data[,1]
  return(strphylo)
}
