#' Generates xml for structured coalescent simulation using MASTER
#'
#' Generates a .xml file to be input to BEAST 2 for structured coalescent simulations using MASTER
#'
#' @param coal_rate vector of coalescent rates
#' @param bit_mig_mat matrix of backward-in-time migration rates
#' @param leaf_data  nx3 matrix or data.frame with first column giving tip labels, second column sample times and third column the initial deme of each sample point
#' @param n_deme number of distinct demes
#' @param xml_path path to output nexus file
#' @param con A connection object or a character string for location to output xml file
#'
#' @return xml file
#'
#' @export

master_xml <- function(coal_rate, bit_mig_mat, leaf_data, n_deme, xml_path, con = stdout()){

  out <- paste("<beast version='2.7' namespace='master:master.model:master.steppers:master.conditions:master.postprocessors:master.outputs'>",
               "<run spec='InheritanceTrajectory' verbosity='0'> \n",
               "<model spec='Model'>",
               paste0("\t <populationType spec='PopulationType' typeName='L' id='L' dim='", n_deme, "'/>"),
               "\t <reactionGroup spec='ReactionGroup' reactionGroupName='Coalescence'>",
               sep = "\n\t")

  for (i in 1 : n_deme){
    out <- paste(out,
                 paste0("<reaction spec='Reaction' rate='", coal_rate[i] / 2, "'>"), #Coalescent rates halved to account for MASTER treating coalescence between (i,j) as distinct from coalescence between (j,i)
                 paste0("\t 2L[", i-1, "]:1 -> L[", i-1, "]:1 "),
                 "</reaction>",
                 sep = "\n\t\t\t")
  }

  out <- paste(out,
               "</reactionGroup> \n",
               "<reactionGroup spec='ReactionGroup' reactionGroupName='Migration'>",
               sep = "\n\t\t")

  for (i in 1 : n_deme){
    for (j in (1:n_deme)[-i]){
      out <- paste(out,
                   paste0("<reaction spec='Reaction' rate='", bit_mig_mat[i,j], "'>"),
                   paste0("\t L[", i-1, "] -> L[", j-1, "]"),
                   "</reaction>",
                   sep = "\n\t\t\t\t")
    }
  }

  out <- paste(out,
               "\t </reactionGroup> \n\n",
               "</model> \n\n",
               "<initialState spec='InitState'>",
               sep = "\n\t")

  leaf_times <- unique(leaf_data[,2])
  n_leaf <- nrow(leaf_data)
  for (time in leaf_times){
    for (deme in 1 : n_deme){
      leaf_set <- (1 : n_leaf)[(leaf_data[,2] == time) & (leaf_data[,3] == deme)]
      if (length(leaf_set) > 0){
        out <- paste(out,
                     paste0("<lineageSeedMultiple spec='MultipleIndividuals' copies='", length(leaf_set), "' time='", time, "'>"),
                     paste0("\t <population spec='Population' type='@L' location='", deme - 1, "'/>"),
                     "</lineageSeedMultiple>",
                     sep = "\n\t\t")
      }
    }
  }

  out <- paste(out,
               "\t </initialState> \n",
               "\t <lineageEndCondition spec='LineageEndCondition' nLineages='1'/> \n",
               paste0("\t\t <output spec='NexusOutput' fileName='", xml_path, "' reverseTime='true'/>"),
               "\t </run>",
               "</beast>",
               sep = "\n")

  writeLines(out, con = con)
}


#' Generates BEAST2 xml for MultiTypeTree
#'
#' Generates a .xml file to run BEAST2 using MultiTypeTree with a fixed dated phylogeny.
#'
#' @param strphylo structured phylo object giving initialisation condition for fixed tree run with MultiType Tree
#' @param coal_rate initial estimate for coalescent rates
#' @param bit_mig_rate initial estimate for backward-in-time migration matrix
#' @param N total number of MCMC iterations (including burn-in)
#' @param thin thinning increment for the MCMC
#' @param migration_history_thin Thinning rate for migration history samples to be saved
#' @param con A connection object or a character string giving the location for the output xml file (stdout() prints to console)
#' @param run_name Name for logger files to be saved as
#' @param priors either 'default' in which lognormal priors are used or 'gamma' in which gamma/inverse gamma priors are used (and prior parameters must be specified)
#' @param proposal_rates Relative weights of proposal operators (migration rates update : effective population size updates : migration history updates )
#'
#' @return output file or file content on screen
#'
#' @export

MultiTypeTree_xml <- function(strphylo, n_deme, coal_rate, bit_mig_mat,
                           N=1e7, thin=1e3, migration_history_thin=thin,
                           con = stdout(), run_name = "$(filebase)",
                           priors = 'default',
                           proposal_rates = c(1,1,1),
                           cr_shape=NULL, cr_rate=NULL,
                           mm_shape=NULL, mm_rate=NULL){

  if (!(priors %in% c('default', 'Default', 'lognormal', 'gamma', 'Gamma'))) stop("Invalid prior distribution selected - use either 'default' for lognormal priors or 'gamma' for gamma priors")
  if ((priors %in% c('gamma', 'Gamma')) && (is.null(cr_shape))) stop('Missing coalescent rate prior parameters')
  if ((priors %in% c('gamma', 'Gamma')) && (is.null(mm_shape))) stop('Missing migration rate prior parameters')

  n_tip <- length(strphylo$tip.label)
  node_ages <- ape::node.depth.edgelength(strphylo)

  cat('<beast version="2.7"',
        'namespace="beast.base.evolution.alignment:',
          '\tbeast.pkgmgmt:',
          '\tbeast.base.core:',
          '\tbeast.base.inference:',
          '\tbeast.pkgmgmt:',
          '\tbeast.base.core:',
          '\tbeast.base.inference.parameter:',
          '\tbeast.base.evolution.tree:',
          '\tbeast.base.evolution.tree.coalescent:',
          '\tbeast.base.util:',
          '\tbeast.base.math:',
          '\tbeast.base.evolution.operator:',
          '\tbeast.base.inference.operator:',
          '\tbeast.base.evolution.sitemodel:',
          '\tbeast.base.evolution.substitutionmodel:',
          '\tbeast.base.evolution.likelihood:',
          '\tbeast.evolution.migrationmodel:',
          '\tbeast.base.inference.distribution:',
          '\tmultitypetree.distributions:',
          '\tmultitypetree.operators:',
          '\tmultitypetree.util">\n',
      file = con,
      sep='\n\t')

  # Genetic sequences for all tips consists of a single '?' base
  cat('\t<!-- Sequence Alignment -->',
        '<alignment spec="beast.base.evolution.alignment.Alignment" id="alignment" dataType="nucleotide">',
          paste0("\t<sequence taxon='", strphylo$tip.label, "' value='?'/>"),
        '</alignment> \n',
      file=con,
      sep = '\n\t',
      append=TRUE)

  cat("\t<!-- Tip sampling demes -->",
      "<typeTraitSet",
      "\tspec='TraitSet'",
      "\tid='typeTraitSet'",
      "\ttraitname='type'",
      "\tvalue='",
      paste0('\t\t', strphylo$tip.label, '=', strphylo$node.deme[1:n_tip], collapse = ',\n\t'),
      "'>",
      "\t<taxa spec='TaxonSet' alignment='@alignment'/>",
      "</typeTraitSet>\n",
      file=con,
      sep = "\n\t",
      append=TRUE)

  cat("\t<!-- Tip sampling times -->",
      "<timeTraitSet",
      "\tspec='TraitSet'",
      "\tid='timeTraitSet'",
      "\ttraitname='date-backward'",
      "\tvalue='",
      paste0('\t\t', strphylo$tip.label, '=', node_ages[1:n_tip], collapse=',\n\t'),
      "'>",
      "\t<taxa spec='TaxonSet' alignment='@alignment'/>",
      "</timeTraitSet>\n",
      file=con,
      sep='\n\t',
      append=TRUE)

  cat("\t<!-- Migration model -->",
      "<migrationModel spec='multitypetree.evolution.tree.SCMigrationModel' id='migModel'>",
      paste0("\t<rateMatrix spec='RealParameter' dimension='", n_deme * (n_deme - 1), "' id=\"rateMatrix\">"),
      paste0("\t\t", paste(bit_mig_mat[-(1 + 0:(n_deme - 1) * (n_deme + 1))], collapse = " ")),
      "\t</rateMatrix>",
      paste0("\t<popSizes spec='RealParameter' dimension=\"", n_deme, "\" id=\"popSizes\">"),
      paste0("\t\t", paste(1/coal_rate, collapse = " ")),
      "\t</popSizes>",
      "\t<typeSet id=\"typeSet\" spec='multitypetree.evolution.tree.TypeSet' typeTraitSet=\"@typeTraitSet\"/>",
      "</migrationModel> \n\n",
      file=con,
      sep='\n\t',
      append=TRUE)

  cat("\t<!-- Parameter priors -->",
      "<input spec='CompoundDistribution' id='parameterPriors'>",
      file=con,
      sep='\n\t',
      append=TRUE)

  if (priors %in% c('default', 'Default', 'lognormal')){ # Default MTT priors - logNormal(0,4)
    cat(# Migration rates
      '\t\t <distribution spec="beast.base.inference.distribution.Prior" x="@rateMatrix">',
      '\t\t <distr spec="LogNormalDistributionModel" M="0.0" S="4.0"/>',
      '\t </distribution>',
      # Coalescent rates
      '\t <distribution spec="beast.base.inference.distribution.Prior" x="@popSizes">',
      '\t\t <distr spec="LogNormalDistributionModel"  M="0.0" S="4.0"/>',
      "\t </distribution>",
      "</input> \n",
      file=con,
      sep = "\n\t",
      append=TRUE)
  } else if (priors %in% c('gamma', 'Gamma')){ #Gamma priors on migration rates and coalescent rates (inverse gamma on effective pop size)
    cat(# Migration rates
      '\t\t <distribution spec="beast.base.inference.distribution.Prior" x="@rateMatrix">',
      paste0('\t\t <distr spec="beast.base.inference.distribution.Gamma" alpha="', mm_shape, '" beta="', mm_rate, '"/>'),
      '\t </distribution>',
      # Coalescent rates
      '\t <distribution spec="beast.base.inference.distribution.Prior" x="@popSizes">',
      paste0('\t\t <distr spec="beast.base.inference.distribution.InverseGamma"  alpha="', cr_shape, '" beta="', cr_rate, '"/>'),
      "\t </distribution>",
      "</input> \n",
      file=con,
      sep = "\n\t",
      append=TRUE)
  }

  cat("\t<!-- Structured coalescent probability density -->",
      "<input spec='StructuredCoalescentTreeDensity' id='treePrior'>",
      "\t <multiTypeTree idref=\"tree\"/>",
      "\t <migrationModel idref=\"migModel\"/>",
      "</input> \n",
      file=con,
      sep='\n\t',
      append=TRUE)

  # Convert strphylo to BEAST2 metacommented Newick
  class(strphylo) <- "phylo"
  treedata <- tidytree::as.treedata(strphylo)
  treedata@data <- tidytree::tibble(type = paste0("\"", strphylo$node.deme - 1, "\""), node = 1:length(strphylo$node.deme))
  newick_tree <- treeio::write.beast.newick(treedata)

  cat('\t<run spec="MCMC" id="mcmc" chainLength="', format(N, scientific=FALSE), '" storeEvery="', format(thin, scientific=FALSE), '">',
      file=con,
      sep='',
      append=TRUE)

  cat("\t\t<!-- Initialise tree from Newick string -->",
      "<init spec='multitypetree.evolution.tree.MultiTypeTreeFromNewick'",
      "\tid='tree' >",
      "\t<![CDATA[",
      paste0('\t\t', newick_tree),
      "\t]]>",
      "\t<typeSet idref='typeSet'/>",
      "\t<typeTrait idref='typeTraitSet'/>",
      "\t<trait idref='timeTraitSet'/>",
      "</init>\n",
      file=con,
      sep='\n\t\t',
      append=TRUE)

  cat('\t\t<!-- Initial state -->',
      '<state>',
      '\t<stateNode idref="tree"/>',
      '\t<stateNode idref="rateMatrix"/>',
      '\t<stateNode idref="popSizes"/>',
      '</state>\n',
      file=con,
      sep='\n\t\t',
      append=TRUE)

  cat("\t\t<!-- Posterior distribution -->",
      "<distribution spec='CompoundDistribution' id='posterior'>",
      "\t<distribution idref='treePrior'/>",
      "\t<distribution idref=\"parameterPriors\"/>",
      "</distribution>\n",
      file=con,
      sep='\n\t\t',
      append=TRUE)

    ##### MTT operators
    cat("\t\t<!-- Migration rates scaler -->",
        "<operator spec='ScaleOperator'",
        "\tid='RateScaler'",
        "\tparameter=\"@rateMatrix\"",
        "\tscaleFactor=\"0.8\"",
        paste0("\tweight=\"", proposal_rates[1], "\"/>\n"),
        "<!-- Effective population sizes scaler -->",
        "<operator spec='ScaleOperator'",
        "\tid='PopSizeScaler'",
        "\tparameter=\"@popSizes\"",
        "\tscaleFactor=\"0.8\"",
        paste0("\tweight=\"", proposal_rates[2], "\"/>\n"),
        "<!-- Node retype operator -->",
        "<operator",
        "\tspec=\"NodeRetype\"",
        "\tid=\"NR\"",
        paste0("\tweight=\"", proposal_rates[3], "\"\n"),
        "\tmultiTypeTree=\"@tree\"",
        "\tmigrationModel=\"@migModel\"/>\n",
        file=con,
        sep = "\n\t\t",
        append=TRUE)

    cat("\t\t<!-- Log file -->",
        "<logger",
        paste0('\tlogEvery="', format(thin, scientific=FALSE), '"'),
        paste0('\tfileName="', run_name, '.log">'),
        '\t<model idref="posterior"/>',
        '\t<log idref="posterior"/>',
        '\t<log idref="treePrior"/>',
        '\t<log spec="TreeHeightLogger" tree="@tree"/>',
        '<log id="migModelLogger" spec="MigrationModelLogger" migrationModel="@migModel" multiTypeTree="@tree"/>',
        "</logger>\n",
        file=con,
        sep='\n\t\t',
        append=TRUE)

    cat('\t\t<!-- Trees file -->',
        "<logger",
        paste0('\tlogEvery="', format(migration_history_thin, scientific=FALSE), '"'),
        paste0('\tfileName="', run_name, '.trees"'),
        '\t mode="tree">',
        '\t <log idref="tree"/>',
        '</logger>\n',
        file=con,
        sep='\n\t\t',
        append=TRUE)

    cat('\t\t<!-- stdout() log -->',
        '<logger logEvery="1000">',
        '\t<model idref="posterior"/>',
        '\t<log idref="posterior"/>',
        '\t<log idref="treePrior"/>',
        '\t<log spec="TreeHeightLogger" tree="@tree"/>',
        '\t<ESS spec="beast.base.inference.util.ESS" name="log" arg="@treePrior"/>',
        '\t <ESS spec="beast.base.inference.util.ESS" name="log" arg="@rateMatrix"/>',
        '\t <ESS spec="beast.base.inference.util.ESS" name="log" arg="@popSizes"/>',
        '</logger>',
        file=con,
        sep = '\n\t\t',
        append=TRUE)

    cat('\t</run>',
        '</beast>',
        file=con,
        sep='\n',
        append=TRUE)
}
