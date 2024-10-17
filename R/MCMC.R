Ewing_2004_mcmc <- function(){
  return(NULL)
}


Vaughan_2015_MultiTypeTree_mcmc <- function(){
  return(NULL)
}

#' StructCoalescent MCMC
#'
#' Runs an MCMC chain using the StructCoalescent method using radius-based subtrees
#'
#' @param N Total number of MCMC iterations to complete
#' @param strphylo Initial structured phylogenetic tree
#' @param coal_rate Initial estimate of coalescent rates
#' @param bit_mig_mat Initial estimate of backward-in-time migration rates matrix
#' @param st_radius Initial subtree radius
#' @param adaptive (logical) Adaptively update subtree radius within run? (default TRUE)
#' @param cr_mode Coalescent rates prior mode
#' @param cr_var Coalescent rates prior variance
#' @param mm_mode Backward-in-time migration rates prior mode
#' @param mm_var Backward-in-time migration rates prior variance
#' @param output_dir Directory to output log files to
#' @param run_name Run name to save log files as
#' @param stdout_log (logical) Output logs to stdout()? Default true
#' @param thin Thinning rate for continuous parameter posterior samples
#' @param save_migration_history Logical value indicating whether tree samples are saved (i.e. whether .trees file is created)
#' @param migration_history_thin Thinning rate for tree samples to be saved
#' @param proposal_rates Relative rates of migration history and continuous parameter updates (migration history : coalescent rates : migration rates)
#' @param adaptation_rate Rate at which adaptive MCMC varies subtree radius (unused if adaptive = FALSE)
#' @param target_accept_rate Target acceptance rate for migration history updates
#'
#' @export

StructCoalescent_mcmc <- function(N,
                                  strphylo, coal_rate, bit_mig_mat,
                                  st_radius, adaptive = TRUE,
                                  cr_mode, cr_var,
                                  mm_mode, mm_var,
                                  output_dir, run_name = 'StructCoalescent', stdout_log = TRUE,
                                  thin = max(N/5e3, 1), save_migration_history = TRUE, migration_history_thin = thin,
                                  proposal_rates=c(1e3, 1, 1),
                                  adaptation_rate = 0.6, target_accept_rate = 0.234){

  ED <- as.ED(strphylo)
  n_deme <- max(ED[,5])

  # Convert priors to rate-shape parameterisation from mode-variance
  cr_rate <- (cr_mode + sqrt(cr_mode^2 + 4 * cr_var))/(2 * cr_var)
  cr_shape <- 1 + cr_mode * cr_rate
  mm_rate <- (mm_mode + sqrt(mm_mode^2 + 4 * mm_var))/(2 * mm_var)
  mm_shape <- 1 + mm_mode * mm_rate

  # Forward-in-time rates and eigen decomposition
  fit_mig_mat <- FitMigMatC(bit_mig_mat, coal_rate)
  fit_rates <- fit_mig_mat
  diag(fit_rates) <- - rowSums(fit_mig_mat)

  eigen_decomp <- eigen(fit_rates)
  eigen_vals <- eigen_decomp$values
  eigen_vecs <- eigen_decomp$vectors
  inverse_vecs <- solve(eigen_vecs)

  # Initial likelihoods and priors
  ED_NI <- NodeIndicesC(ED)
  ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)
  mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
  cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))


  # Write iteration number, likelihood, posterior and subtree radius to stdout()
  cat('Sample',
      'Likelihood',
      'Posterior',
      'Subtree radius\n',
      sep = '\t')

  cat(0, # Sample
      sprintf('%.03f', ED_SC), # Likelihood
      sprintf('%.03f', ED_SC + mm_prior + cr_prior), # posterior
      sprintf('%.03f\n', st_radius), #subtree radius (being adapted!)
      sep = '\t')

  # Set up .freq file to store move acceptance frequencies
  freq_file <- file.path(output_dir,
                         paste0(run_name, '.freq'))
  freq <- matrix(0, 2, 3,
                 dimnames = list(c("#accept", "#reject"), c("DTA", "CR", "MM")))

  # Set up .log file to store posterior continuous parameters
  log_file <- file.path(output_dir,
                        paste0(run_name, '.log'))
  cat("sample",
      "likelihood",
      "posterior",
      paste("coal_rate_", 1:n_deme, sep = "", collapse =","),
      paste("backward_migration_rate_", as.vector(outer(1:n_deme, 1:n_deme, paste, sep = "_")[-(1 + 0:(n_deme - 1) * (n_deme + 1))]), sep = ""),
      "st_radius",
      file = log_file, sep =",")

  cat(paste0("\n", 0), #sample
      ED_SC, # likelihood
      ED_SC + mm_prior + cr_prior, # posterior
      coal_rate, #coal_rate
      as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #backward_migration_rate
      st_radius, #Subtree radius
      file = log_file, sep =",", append = TRUE)

  # Set up .trees file to store posterior sampled trees
  if (save_migration_history){
    tree_file <- file.path(output_dir,
                           paste0(run_name, '.trees'))
    tree_data <- as.treedata(ED)
    header <- capture.output(treeio::write.beast(tree_data, file=stdout(), translate=TRUE, tree.name='STATE_0')) #Generate full .trees file for initial tree - need all except final "END;"
    header[2] <- paste("[R-package StructCoalescent, ", date(), "]\n\n", sep = "") #Update package line of .trees file to scoal
    cat(header[-length(header)], file = tree_file, sep = "\n") #Save file with updated package line, omitting "END;" on final line
  }

  for (x in 1 : N){
    move_id <- sample(1:3, 1, prob = proposal_rates)

    if (move_id == 1){
      subtree <- st_centre_dist(ED, st_radius, ED_NI)
      proposal <- local_DTA_subtree_proposal(subtree$ED, subtree$st_labels, fit_rates,
                                             eigen_decomp = eigen_decomp, inverse_vecs = inverse_vecs)
      prop <- proposal$proposal

      if ((nrow(ED) == nrow(prop)) &&(all(na.omit(as.vector(ED == prop))))){ # Early acceptance if prop == EED
        freq[1, move_id] <- freq[1, move_id] + 1
        log_AR <- 0
      } else {
        prop_NI <- NodeIndicesC(prop)
        prop_SC <- SC_like_C(prop, coal_rate, bit_mig_mat, prop_NI)

        log_AR <- min(0, prop_SC - ED_SC +
                        local_DTA_likelihood(subtree$st_labels, coal_rate, bit_mig_mat)$log.likelihood -
                        proposal$prop_prob)

        if (log(runif(1)) < log_AR){
          # Accept
          freq[1, move_id] <- freq[1, move_id] + 1

          ED <- prop
          ED_NI <- prop_NI
          ED_SC <- prop_SC
        }
      }

      if (adaptive){
        # Update proposal radius
        st_radius <- exp(log(st_radius) + x^(-adaptation_rate) * (exp(log_AR) - target_accept_rate))
      }
    } else {
      if (move_id == 2){
        coal_rate <- ED_cr_gibbs(ED, n_deme, ED_NI, cr_shape, cr_rate)
        cr_prior <- sum(dgamma(coal_rate, cr_shape, cr_rate, log = TRUE))
      } else if (move_id == 3){
        bit_mig_mat <- ED_bmm_gibbs(ED, n_deme, ED_NI, mm_shape, mm_rate)
        mm_prior <- sum(dgamma(bit_mig_mat, mm_shape, mm_rate, log = TRUE)[-(1 + 0:(n_deme - 1) * (n_deme + 1))])
      }

      freq[1, move_id] <- freq[1, move_id] + 1 #Gibbs move always accepted

      # Update forward-in-time rates, pre-computed likelihoods and eigendecomposition
      ED_SC <- SC_like_C(ED, coal_rate, bit_mig_mat, ED_NI)

      fit_mig_mat <- FitMigMatC(bit_mig_mat, coal_rate)
      fit_rates <- fit_mig_mat
      diag(fit_rates) <- - rowSums(fit_mig_mat)

      eigen_decomp <- eigen(fit_rates)
      eigen_vals <- eigen_decomp$values
      eigen_vecs <- eigen_decomp$vectors
      inverse_vecs <- solve(eigen_vecs)
    }

    freq[2, move_id] <- freq[2, move_id] + 1 # Increment proposal quantity for move_id

    if (x %% thin == 0){
      # Write iteration, likelihood, posterior and subtree radius to stdout()
      cat(x, # Sample
          sprintf('%.03f', ED_SC), # Likelihood
          sprintf('%.03f', ED_SC + mm_prior + cr_prior), # posterior
          sprintf('%.03f\n', st_radius), #subtree radius (being adapted!)
          sep = '\t')

      # Write continuous parameters to .log file
      cat(paste0("\n", x), #sample
          ED_SC, # likelihood
          ED_SC + mm_prior + cr_prior, # posterior
          coal_rate, #coal_rate
          as.vector(bit_mig_mat)[-(1 + 0:(n_deme - 1) * (n_deme + 1))], #Mig mat
          st_radius, #subtree radius (being adapted!)
          file = log_file, sep =",", append = TRUE)

      # Update .freq file (Overwrites existing file entirely)
      write.table(freq,
                  file = freq_file,
                  row.names = c('#ACCEPT', '#TOTAL'),
                  col.names = TRUE)
    }

    if (save_migration_history && (x %% migration_history_thin == 0)){
      # Write current tree to .trees file
      tree_data <- as.treedata(ED)
      cat("\tTREE STATE_", x, " = ",
          treeio::write.beast.newick(tree_data), "\n",
          file = tree_file, append = TRUE, sep = "")
    }
  }
  cat("END;\n", file = tree_file, append = TRUE, sep = "")
}
