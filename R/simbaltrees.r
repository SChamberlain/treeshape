#' Simulate a set of balanced and unbalanced trees. 
#' 
#' Could use this in asking questions aobut how phylogenetic tree balance 
#'     influences ____ (madlib it). 
#'     
#' The Ornstein-Uhlenbeck code is commented out within the function right now
#' 		as it was replaced pretty more or less with the early burst (EB) model.
#' 		
#' @import ape plyr apTreeshape bipartite reshape2 geiger
#' @importFrom phytools fastBM
#' @importFrom phytools pbtree
#' @importFrom phytools nodeHeights
#' @param tips Number of species to simulate in each tree - will be same for 
#' 		all trees.
#' @param metric Methods to use to generate trees, one of "colless", "beta", or  
#' 		gamma (see details). Defaults to "colless".
#' @param numtrees Number of trees to produce. Defaults to 10 trees.
#' @param cutlow Value at which to filter trees on the low (e.g., unbalanced) 
#' 		side of the metric.
#' @param cuthigh Value at which to filter trees on the high (e.g., balanced) 
#' 		side of the metric.
#' @param alpha Alpha value for the Ornstein-Uhlenbeck model of trait evolution.  
#' 		From ape documentation: "a numeric vector giving the strength of the 
#' 		selective constraint for each branch (can be a single value)."
#' @param sigma Sigma is the single value of the standard-deviation of the 
#' 		random component for each branch (can be a single value).
#' @param netmets Network structure metrics to calculate - only use those that 
#' 		calculate single values for each matrix.
#' @param alpha_eb Alpha parameter for the Early Burst trait simulation model. 
#' 		From the geiger documentation: "Exponent of the relationship between 
#' 		rate and time in the exponentialchange model".
#' @param sigma_eb Sigma parameter for the Early Burst trait simulation model. 
#' 		From the ape documentation: "Sigma is the single value of the 
#' 		standard-deviation of the random component for each branch (can be a 
#' 		single value).".  This sigma means the same as for the OU model, but this
#' 		allows to specify it separately.
#' @param rval Value for the ratio model for constructing interaction network 
#' 		matrices.
#' @param cval Value for the complimentarity model for constructing interation
#' 		network matrices. 
#' @return A data.frame of network structure metrics for balanced and unbalanced 
#' 		trees.
#' @examples \dontrun{
#' netmets <- c("connectance", "links per species", "nestedness", "web asymmetry")
#' temp <- simbaltrees(tips=10, metric="colless", numtrees=10, cutlow=-0.5, cuthigh=0.5, alpha_eb=-0.8, sigma_eb=3, rval=1.5, cval=2, netmets=netmets)
#' }
#' @export
simbaltrees <- function(tips = 10, metric, numtrees, cutlow, cuthigh, alpha, 
	sigma, alpha_eb, sigma_eb, rval, cval,
	netmets = c("connectance", "links per species", "nestedness", "web asymmetry")) 
{  
  trees_colless_plants <- simbal(t=tips, metric=metric, n=numtrees, 
  															 cutlow = cutlow, cuthigh = cuthigh)
  trees_colless_plants_bal <- trees_colless_plants$bal # get the balanced trees
  trees_colless_plants_unbal <- trees_colless_plants$unbal # get the unbalanced trees
  
  ### animal trees
  trees_colless_anim <- simbal(t=tips, metric=metric, n=numtrees, 
  														 cutlow = cutlow, cuthigh = cuthigh)
  trees_colless_anim_bal <- trees_colless_anim$bal # get the balanced trees
  trees_colless_anim_unbal <- trees_colless_anim$unbal # get the unbalanced trees
  
  ################## Simulate traits on each tree
  # Plants
  ## Brownian motion traits
  t1_col_plants_bal_bm <- lapply(trees_colless_plants_bal, fastBM, a = 10, 
  															 bounds=c(0,100))
  t1_col_plants_unbal_bm <- lapply(trees_colless_plants_unbal, fastBM, a = 10, 
  																 bounds=c(0,100))
  
  ## Random traits - randomize traits across tips - shuffles names, but doesn't move position in tree of course
  t1_col_plants_bal_rand <- lapply(t1_col_plants_bal_bm, randtips)
  t1_col_plants_unbal_rand <- lapply(t1_col_plants_unbal_bm, randtips)
  
  ## Orntsein-Uhlenbeck model for conservation of trait evolution - 1 optima
#   t1_col_plants_bal_ou <- lapply(trees_colless_plants_bal, rTraitCont, 
#   															 model = "OU", sigma = sigma, alpha=alpha, theta=1)
#   t1_col_plants_unbal_ou <- lapply(trees_colless_plants_unbal, rTraitCont, 
#   																 model = "OU", sigma = sigma, alpha=alpha, theta=1)
  
  ## Early Burst model
  t1_col_plants_bal_eb <- lapply(trees_colless_plants_bal, function(x) 
  	rTraitCont(exponentialchangeTree(x, a=alpha_eb), sigma = sigma_eb))
  t1_col_plants_unbal_eb <- lapply(trees_colless_plants_unbal, function(x) 
  	rTraitCont(exponentialchangeTree(x, a=alpha_eb), sigma = sigma_eb))
  
  # Animals
  ## Brownian motion traits
  t1_col_anim_bal_bm <- llply(trees_colless_anim_bal, fastBM, a = 10, 
  														bounds=c(0,100))
  t1_col_anim_unbal_bm <- llply(trees_colless_anim_unbal, fastBM, a = 10, 
  															bounds=c(0,100))
  
  ## Random traits - randomize traits across tips - shuffles names, but doesn't move position in tree of course
  t1_col_anim_bal_rand <- lapply(t1_col_anim_bal_bm, randtips)
  t1_col_anim_unbal_rand <- lapply(t1_col_anim_unbal_bm, randtips)
  
  ## Orntsein-Uhlenbeck model for conservation of trait evolution - 1 optima
#   t1_col_anim_bal_ou <- lapply(trees_colless_anim_bal, rTraitCont, 
#   														 model = "OU", sigma = sigma, alpha=alpha, theta=1)
#   t1_col_anim_unbal_ou <- lapply(trees_colless_anim_unbal, rTraitCont, 
#   															 model = "OU", sigma = sigma, alpha=alpha, theta=1)
  
  ## Early Burst model
  t1_col_anim_bal_eb <- lapply(trees_colless_anim_bal, function(x) 
  	rTraitCont(exponentialchangeTree(x, a=alpha_eb), sigma = sigma_eb))
  t1_col_anim_unbal_eb <- lapply(trees_colless_anim_unbal, function(x) 
  	rTraitCont(exponentialchangeTree(x, a=alpha_eb), sigma = sigma_eb))
  
  
  ################## Measure aspects of traits on trees
  # Plants
  t_p_bal_bm <- traitsig(t1_col_plants_bal_bm, trees_colless_plants_bal)
  t_p_bal_rand <- traitsig(t1_col_plants_bal_rand, trees_colless_plants_bal)
#   t_p_bal_ou <- traitsig(t1_col_plants_bal_ou, trees_colless_plants_bal)
  t_p_bal_eb <- traitsig(t1_col_plants_bal_eb, trees_colless_plants_bal)
  t_p_unbal_bm <- traitsig(t1_col_plants_unbal_bm, trees_colless_plants_unbal)
  t_p_unbal_rand <- traitsig(t1_col_plants_unbal_rand, trees_colless_plants_unbal)
#   t_p_unbal_ou <- traitsig(t1_col_plants_unbal_ou, trees_colless_plants_unbal)
  t_p_unbal_eb <- traitsig(t1_col_plants_unbal_eb, trees_colless_plants_unbal)
  
  # Animals
  t_a_bal_bm <- traitsig(t1_col_anim_bal_bm, trees_colless_anim_bal)
  t_a_bal_rand <- traitsig(t1_col_anim_bal_rand, trees_colless_anim_bal)
#   t_a_bal_ou <- traitsig(t1_col_anim_bal_ou, trees_colless_anim_bal)
  t_a_bal_eb <- traitsig(t1_col_anim_bal_eb, trees_colless_anim_bal)
  t_a_unbal_bm <- traitsig(t1_col_anim_unbal_bm, trees_colless_anim_unbal)
  t_a_unbal_rand <- traitsig(t1_col_anim_unbal_rand, trees_colless_anim_unbal)
#   t_a_unbal_ou <- traitsig(t1_col_anim_unbal_ou, trees_colless_anim_unbal)
  t_a_unbal_eb <- traitsig(t1_col_anim_unbal_eb, trees_colless_anim_unbal)
  
  # Make data.frame of results
  traits_list <- list(t_p_bal_bm, t_p_bal_rand, t_p_bal_eb, 
  										t_p_unbal_bm, t_p_unbal_rand, t_p_unbal_eb,
                      t_a_bal_bm, t_a_bal_rand, t_a_bal_eb, 
  										t_a_unbal_bm, t_a_unbal_rand, t_a_unbal_eb)
  names(traits_list) <- c("t_p_bal_bm", "t_p_bal_rand", "t_p_bal_eb", 
  												"t_p_unbal_bm", "t_p_unbal_rand",
  												"t_p_unbal_eb", "t_a_bal_bm", "t_a_bal_rand", 
  												"t_a_bal_eb", "t_a_unbal_bm", "t_a_unbal_rand", 
  												 "t_a_unbal_eb")
  traits_df <- ldply(traits_list)
  traits_df$numsp <- rep(tips, nrow(traits_df))
  
  ################## Get plant-animal phylogenetic tree pairs
  # Balanced trees
  tree_pairs_bal <- list(trees_colless_plants_bal, trees_colless_anim_bal)
  
  # Unbalanced trees
  tree_pairs_unbal <- list(trees_colless_plants_unbal, trees_colless_anim_unbal)
  
  # Pairs of traits from balanced trees
  ## Brownian motion traits
  all_t1_bal_bm <- list(t1_col_plants_bal_bm, t1_col_anim_bal_bm)
 
  ## Random traits
  all_t1_bal_rand <- list(t1_col_plants_bal_rand, t1_col_anim_bal_rand)
 
  ## OU traits
#   all_t1_bal_ou <- list(t1_col_plants_bal_ou, t1_col_anim_bal_ou)
  
  ## EB traits
  all_t1_bal_eb <- list(t1_col_plants_bal_eb, t1_col_anim_bal_eb)

  # Pairs of traits form unbalanced trees
  ## Brownian motion traits
  all_t1_unbal_bm <- list(t1_col_plants_unbal_bm, t1_col_anim_unbal_bm)
 
  ## Random traits
  all_t1_unbal_rand <- list(t1_col_plants_unbal_rand, t1_col_anim_unbal_rand)
 
  ## OU traits
#   all_t1_unbal_ou <- list(t1_col_plants_unbal_ou, t1_col_anim_unbal_ou)
  
  ## EB traits
  all_t1_unbal_eb <- list(t1_col_plants_unbal_eb, t1_col_anim_unbal_eb)
  
  ################## Simulate networks on each tree
  # Simulate completely random networks, regardless of traits, etc.
  mats_rand_bal <- sim_rand_nets(tree_pairs_bal)
  mats_rand_unbal <- sim_rand_nets(tree_pairs_unbal)
  
  # Simulate networks, with interactions propoprtional to trait matching
  #   listoftraitvecs <- all_t1_bal_ou
  # ratio
  ## BM
  mats_traits1_bal_bm_ratio <- sim_traits_nets(all_t1_bal_bm, method = "r", value = rval)
  mats_traits1_unbal_bm_ratio <- sim_traits_nets(all_t1_unbal_bm, method = "r", value = rval)
  
  ## Random
  mats_traits1_bal_rand_ratio <- sim_traits_nets(all_t1_bal_rand, method = "r", value = rval)
  mats_traits1_unbal_rand_ratio <- sim_traits_nets(all_t1_unbal_rand, method = "r", value = rval)
  
  ## OU
#   mats_traits1_bal_ou_ratio <- sim_traits_nets(all_t1_bal_ou, method = "r", value = rval)
#   mats_traits1_unbal_ou_ratio <- sim_traits_nets(all_t1_unbal_ou, method = "r", value = rval)
  
  ## EB
  mats_traits1_bal_eb_ratio <- sim_traits_nets(all_t1_bal_eb, method = "r", value = rval)
  mats_traits1_unbal_eb_ratio <- sim_traits_nets(all_t1_unbal_eb, method = "r", value = rval)
  
  # complementarity	
  ## BM
  mats_traits1_bal_bm_comp <- sim_traits_nets(all_t1_bal_bm, method = "c", value = cval)
  mats_traits1_unbal_bm_comp <- sim_traits_nets(all_t1_unbal_bm, method = "c", value = cval)
 
  ## Random
  mats_traits1_bal_rand_comp <- sim_traits_nets(all_t1_bal_rand, method = "c", value = cval)
  mats_traits1_unbal_rand_comp <- sim_traits_nets(all_t1_unbal_rand, method = "c", value = cval)

  ## OU
#   mats_traits1_bal_ou_comp <- sim_traits_nets(all_t1_bal_ou, method = "c", value = cval)
#   mats_traits1_unbal_ou_comp <- sim_traits_nets(all_t1_unbal_ou, method = "c", value = cval)
  
  ## EB
  mats_traits1_bal_eb_comp <- sim_traits_nets(all_t1_bal_eb, method = "c", value = cval)
  mats_traits1_unbal_eb_comp <- sim_traits_nets(all_t1_unbal_eb, method = "c", value = cval)
  
  # barrier (no need to give value parameter)
  ## BM
  mats_traits1_bal_bm_barr <- sim_traits_nets(all_t1_bal_bm, method = "b")
  mats_traits1_unbal_bm_barr <- sim_traits_nets(all_t1_unbal_bm, method = "b")
 
  ## Random
  mats_traits1_bal_rand_barr <- sim_traits_nets(all_t1_bal_rand, method = "b")
  mats_traits1_unbal_rand_barr <- sim_traits_nets(all_t1_unbal_rand, method = "b")

  ## OU
#   mats_traits1_bal_ou_barr <- sim_traits_nets(all_t1_bal_ou, method = "b")
#   mats_traits1_unbal_ou_barr <- sim_traits_nets(all_t1_unbal_ou, method = "b")
  
  ## EB
  mats_traits1_bal_eb_barr <- sim_traits_nets(all_t1_bal_eb, method = "b")
  mats_traits1_unbal_eb_barr <- sim_traits_nets(all_t1_unbal_eb, method = "b")
  
  ################## Calculate network metrics on matrices
  df_rand <- getnetmets(mats_rand_bal, mats_rand_unbal, netmets=netmets) # random networks
  
  ## BM
  df_traits1_bm_ratio <- getnetmets(mats_traits1_bal_bm_ratio, mats_traits1_unbal_bm_ratio, netmets=netmets) # ratio traits
  df_traits1_bm_comp <- getnetmets(mats_traits1_bal_bm_comp, mats_traits1_unbal_bm_comp, netmets=netmets) # complementarity traits
  df_traits1_bm_barr <- getnetmets(mats_traits1_bal_bm_barr, mats_traits1_unbal_bm_barr, netmets=netmets) # barrier traits
  
  ## Random
  df_traits1_rand_ratio <- getnetmets(mats_traits1_bal_rand_ratio, mats_traits1_unbal_rand_ratio, netmets=netmets) # ratio traits
  df_traits1_rand_comp <- getnetmets(mats_traits1_bal_rand_comp, mats_traits1_unbal_rand_comp, netmets=netmets) # complementarity traits
  df_traits1_rand_barr <- getnetmets(mats_traits1_bal_rand_barr, mats_traits1_unbal_rand_barr, netmets=netmets) # barrier traits
  
  ## OU
#   df_traits1_ou_ratio <- getnetmets(mats_traits1_bal_ou_ratio, mats_traits1_unbal_ou_ratio, netmets=netmets) # ratio traits
#   df_traits1_ou_comp <- getnetmets(mats_traits1_bal_ou_comp, mats_traits1_unbal_ou_comp, netmets=netmets) # complementarity traits
#   df_traits1_ou_barr <- getnetmets(mats_traits1_bal_ou_barr, mats_traits1_unbal_ou_barr, netmets=netmets) # barrier traits
  
  ## EB
  df_traits1_eb_ratio <- getnetmets(mats_traits1_bal_eb_ratio, mats_traits1_unbal_eb_ratio, netmets=netmets) # ratio traits
  df_traits1_eb_comp <- getnetmets(mats_traits1_bal_eb_comp, mats_traits1_unbal_eb_comp, netmets=netmets) # complementarity traits
  df_traits1_eb_barr <- getnetmets(mats_traits1_bal_eb_barr, mats_traits1_unbal_eb_barr, netmets=netmets) # barrier traits
  
  alldat <- rbind(df_rand, 
                  df_traits1_bm_ratio, df_traits1_bm_comp, df_traits1_bm_barr,
                  df_traits1_rand_ratio, df_traits1_rand_comp, df_traits1_rand_barr,
  								df_traits1_eb_ratio, df_traits1_eb_comp, df_traits1_eb_barr)
  alldat$traitevol <- c( rep("Random",nrow(df_rand)), 
                         rep("BM",nrow(df_traits1_bm_ratio)),
                         rep("BM",nrow(df_traits1_bm_comp)),
                         rep("BM",nrow(df_traits1_bm_barr)),
                         rep("Random",nrow(df_traits1_rand_ratio)),
                         rep("Random",nrow(df_traits1_rand_comp)),
                         rep("Random",nrow(df_traits1_rand_barr)),
  											 rep("EB",nrow(df_traits1_eb_ratio)),
  											 rep("EB",nrow(df_traits1_eb_comp)),
  											 rep("EB",nrow(df_traits1_eb_barr))
  											)
  alldat$netmethod <- c( rep("Random",nrow(df_rand)), 
                         rep("Ratio",nrow(df_traits1_bm_ratio)),
                         rep("Complementarity",nrow(df_traits1_bm_comp)),
                         rep("Barrier",nrow(df_traits1_bm_barr)),
                         rep("Ratio",nrow(df_traits1_rand_ratio)),
                         rep("Complementarity",nrow(df_traits1_rand_comp)),
                         rep("Barrier",nrow(df_traits1_rand_barr)),
  											 rep("Ratio",nrow(df_traits1_eb_ratio)),
  											 rep("Complementarity",nrow(df_traits1_eb_comp)),
  											 rep("Barrier",nrow(df_traits1_eb_barr))
  											)
  alldat$numsp <- rep(tips, nrow(alldat))
  list(traits_df, alldat) 
}