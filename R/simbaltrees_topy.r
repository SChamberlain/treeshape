#' Simulate a set of balanced and unbalanced trees. 
#' 
#' Could use this in asking questions aobut how phylogenetic tree balance 
#'     influences ____ (madlib it). 
#' 		
#' @import ape plyr apTreeshape bipartite reshape2 geiger
#' @param tips_p Number of plant species to simulate in each tree - will be same for 
#' 		all trees. The number of animal species will be calculated based on the value
#' 		of the asymm parameter. 
#' @param metric Methods to use to generate trees, one of "colless", "beta", or  
#' 		gamma (see details). Defaults to "colless".
#' @param numtrees Number of trees to produce. Defaults to 10 trees.
#' @param cutlow Value at which to filter trees on the low (e.g., unbalanced) 
#' 		side of the metric.
#' @param cuthigh Value at which to filter trees on the high (e.g., balanced) 
#' 		side of the metric.
#' @param a If model = "bm", a value for ancestral state at the root node.
#' @param bounds If model = "bm", a vector with the lower and upper bounds 
#' 		(respectively) for bounded Brownian simulation - by default simulation is unbounded.
#' @param alpha If model = "ou", a value for the Ornstein-Uhlenbeck model of trait evolution.  
#' 		From ape documentation: "a numeric vector giving the strength of the 
#' 		selective constraint for each branch (can be a single value)."
#' @param sigma If model = "ou", is the single value of the standard-deviation of the 
#' 		random component for each branch (can be a single value).
#' @param theta If model = "ou", a numeric vector giving the optimum for each branch (can be 
#' 		a single value)
#' @param alpha_eb If model = "eb", is the exponent of the relationship between 
#' 		rate and time in the exponentialchange model.
#' @param sigma_eb If model = "eb", the single value of the 
#' 		standard-deviation of the random component for each branch (can be a 
#' 		single value).  This sigma means the same as for the OU model, but this
#' 		allows to specify it separately.
#' @param rval Value for the ratio model for constructing interaction network 
#' 		matrices.
#' @param cval Value for the complimentarity model for constructing interation
#' 		network matrices. 
#' @param asymm An asymmetry value in terms of ratio of animals:plants instead of 
#' 		the traditional definition like (Na-Np/Na+Np), where Na is number of animal
#' 		species, and Np is number of plant species. For example, asymm=2 would mean
#' 		twice as many animal as plant species.
#' @param cores Number of cores to use in mcmapply in simulating networks form traits.
#' @param dumpmatrices If FALSE (default) matrices are not spit out to file, but
#' 		if TRUE, then matrices are written to file in separate folders for each (logical)
#' @param matdir Directory to output matrices to. Ignored if dumpmatrices=FALSE.
#' @param netmets Network structure metrics to calculate - only use those that 
#' 		calculate single values for each matrix.
#' @return A data.frame of network structure metrics for balanced and unbalanced 
#' 		trees.
#' @examples \dontrun{
#' temp <- simbaltrees_topy(tips_p=15, metric="colless", numtrees=5, cutlow=-0.5, cuthigh=0.5, a=10, bounds=c(0,100), alpha=1, sigma=1, alpha_eb=-0.8, sigma_eb=3, cval=0.5, asymm=2, matdir="~/newfiles2", modeltorun="complementarity")
#' head(temp) # traits data.frame
#' }
#' @export
simbaltrees_topy <- function(tips_p = 10, metric, numtrees, cutlow, cuthigh, a, 
		bounds, alpha, sigma, theta, alpha_eb, sigma_eb, rval = NULL, cval = NULL, asymm=1, cores = 2,
		dumpmatrices=TRUE, matdir = "~", modeltorun="complementarity") 
{
	### Calculate number of species for plants and animals
	tips_a <- round(tips_p * asymm, 0)
	
	message("Simulating with tips_p=", tips_p, " and tips_a=", tips_a, "...")
	message("...simulating trees...")
	### plant trees
	trees_colless_plants <- simbal(t=tips_p, metric=metric, n=numtrees, 
																 cutlow = cutlow, cuthigh = cuthigh)
	trees_colless_plants_bal <- trees_colless_plants$bal # get the balanced trees
	trees_colless_plants_unbal <- trees_colless_plants$unbal # get the unbalanced trees
	
	### animal trees
	trees_colless_anim <- simbal(t=tips_a, metric=metric, n=numtrees, 
															 cutlow = cutlow, cuthigh = cuthigh)
	trees_colless_anim_bal <- trees_colless_anim$bal # get the balanced trees
	trees_colless_anim_unbal <- trees_colless_anim$unbal # get the unbalanced trees
	
	################## Simulate traits on each tree
	message("...simulating traits on trees...")
	# Plants
	## Brownian motion traits
	t1_col_plants_bal_bm <- sim_traits_ontrees(trees_colless_plants_bal, "bm", a=a, bounds=bounds)
	t1_col_plants_unbal_bm <- sim_traits_ontrees(trees_colless_plants_unbal, "bm", a=a, bounds=bounds)
	
	## Orntsein-Uhlenbeck model for conservation of trait evolution - 1 optima
	t1_col_plants_bal_ou <- sim_traits_ontrees(trees_colless_plants_bal, "ou", alpha=alpha, sigma=sigma, theta=1)
	t1_col_plants_unbal_ou <- sim_traits_ontrees(trees_colless_plants_unbal, "ou", alpha=alpha, sigma=sigma, theta=1)
	
	## Early Burst model
	t1_col_plants_bal_eb <- sim_traits_ontrees(trees_colless_plants_bal, "eb", alpha_eb=alpha_eb, sigma_eb=sigma_eb)
	t1_col_plants_unbal_eb <- sim_traits_ontrees(trees_colless_plants_unbal, "eb", alpha_eb=alpha_eb, sigma_eb=sigma_eb)
	
	# Animals
	## Brownian motion traits
	t1_col_anim_bal_bm <- sim_traits_ontrees(trees_colless_anim_bal, "bm", a=a, bounds=bounds)
	t1_col_anim_unbal_bm <- sim_traits_ontrees(trees_colless_anim_unbal, "bm", a=a, bounds=bounds)
	
	## Orntsein-Uhlenbeck model for conservation of trait evolution - 1 optima
	t1_col_anim_bal_ou <- sim_traits_ontrees(trees_colless_anim_bal, "ou", alpha=alpha, sigma=sigma, theta=1)
	t1_col_anim_unbal_ou <- sim_traits_ontrees(trees_colless_anim_unbal, "ou", alpha=alpha, sigma=sigma, theta=1)
	
	## Early Burst model
	t1_col_anim_bal_eb <- sim_traits_ontrees(trees_colless_anim_bal, "eb", alpha_eb=alpha_eb, sigma_eb=sigma_eb)
	t1_col_anim_unbal_eb <- sim_traits_ontrees(trees_colless_anim_unbal, "eb", alpha_eb=alpha_eb, sigma_eb=sigma_eb)
	
	################## Measure aspects of traits on trees
	message("...measuring trait characteristics...")
	# Plants
	t_p_bal_bm <- traitsig(t1_col_plants_bal_bm, trees_colless_plants_bal)
	t_p_bal_ou <- traitsig(t1_col_plants_bal_ou, trees_colless_plants_bal)
	t_p_bal_eb <- traitsig(t1_col_plants_bal_eb, trees_colless_plants_bal)
	t_p_unbal_bm <- traitsig(t1_col_plants_unbal_bm, trees_colless_plants_unbal)
	t_p_unbal_ou <- traitsig(t1_col_plants_unbal_ou, trees_colless_plants_unbal)
	t_p_unbal_eb <- traitsig(t1_col_plants_unbal_eb, trees_colless_plants_unbal)
	
	# Animals
	t_a_bal_bm <- traitsig(t1_col_anim_bal_bm, trees_colless_anim_bal)
	t_a_bal_ou <- traitsig(t1_col_anim_bal_ou, trees_colless_anim_bal)
	t_a_bal_eb <- traitsig(t1_col_anim_bal_eb, trees_colless_anim_bal)
	t_a_unbal_bm <- traitsig(t1_col_anim_unbal_bm, trees_colless_anim_unbal)
	t_a_unbal_ou <- traitsig(t1_col_anim_unbal_ou, trees_colless_anim_unbal)
	t_a_unbal_eb <- traitsig(t1_col_anim_unbal_eb, trees_colless_anim_unbal)
	
	# Make data.frame of results
	traits_list <- list(t_p_bal_bm, t_p_bal_ou, t_p_bal_eb, 
											t_p_unbal_bm, t_p_unbal_ou, t_p_unbal_eb,
											t_a_bal_bm, t_a_bal_ou, t_a_bal_eb, 
											t_a_unbal_bm, t_a_unbal_ou, t_a_unbal_eb)
	names(traits_list) <- c("t_p_bal_bm", "t_p_bal_ou", "t_p_bal_eb", 
													"t_p_unbal_bm", "t_p_unbal_ou",
													"t_p_unbal_eb", "t_a_bal_bm", "t_a_bal_ou", 
													"t_a_bal_eb", "t_a_unbal_bm", "t_a_unbal_ou", 
													"t_a_unbal_eb")
	traits_df <- ldply(traits_list)
	traits_df$numsp_p <- rep(tips_p, nrow(traits_df))
	traits_df$numsp_a <- rep(tips_a, nrow(traits_df))
	traits_df$numsp_all <- rep(sum(tips_a,tips_p), nrow(traits_df))
	
	################## Get plant-animal phylogenetic tree pairs
# 	# Balanced trees
# 	tree_pairs_bal <- list(trees_colless_plants_bal, trees_colless_anim_bal)
# 	
# 	# Unbalanced trees
# 	tree_pairs_unbal <- list(trees_colless_plants_unbal, trees_colless_anim_unbal)
	
	# Pairs of traits from balanced trees
	## Brownian motion traits
	all_t1_bal_bm <- list(t1_col_plants_bal_bm, t1_col_anim_bal_bm)
	
	## OU traits
	all_t1_bal_ou <- list(t1_col_plants_bal_ou, t1_col_anim_bal_ou)
	
	## EB traits
	all_t1_bal_eb <- list(t1_col_plants_bal_eb, t1_col_anim_bal_eb)
	
	# Pairs of traits form unbalanced trees
	## Brownian motion traits
	all_t1_unbal_bm <- list(t1_col_plants_unbal_bm, t1_col_anim_unbal_bm)
	
	## OU traits
	all_t1_unbal_ou <- list(t1_col_plants_unbal_ou, t1_col_anim_unbal_ou)
	
	## EB traits
	all_t1_unbal_eb <- list(t1_col_plants_unbal_eb, t1_col_anim_unbal_eb)
	
	################## Simulate networks on each tree
	message("...simulating networks using traits and writing to files...")
	# Simulate completely random networks, regardless of traits, etc.
# 	mats_rand_bal <- sim_rand_nets(tree_pairs_bal)
# 	mats_rand_unbal <- sim_rand_nets(tree_pairs_unbal)
	
	if(modeltorun=="complementarity"){
		# complementarity	
		## BM
		sim_traits_nets_par(all_t1_bal_bm, method = "c", value = cval, type="bal", traitm="bm", matdir=matdir)
		sim_traits_nets_par(all_t1_unbal_bm, method = "c", value = cval, type="unbal", traitm="bm", matdir=matdir)
		
		## OU
		sim_traits_nets_par(all_t1_bal_ou, method = "c", value = cval, type="bal", traitm="ou", matdir=matdir)
		sim_traits_nets_par(all_t1_unbal_ou, method = "c", value = cval, type="unbal", traitm="ou", matdir=matdir)
		
		## EB
		sim_traits_nets_par(all_t1_bal_eb, method = "c", value = cval, type="bal", traitm="eb", matdir=matdir)
		sim_traits_nets_par(all_t1_unbal_eb, method = "c", value = cval, type="unbal", traitm="eb", matdir=matdir)
	} else
		if(modeltorun=="barrier"){
			# barrier (no need to give value parameter)
			## BM
			sim_traits_nets_par(all_t1_bal_bm, method = "b", type="bal", traitm="bm", matdir=matdir)
			sim_traits_nets_par(all_t1_unbal_bm, method = "b", type="unbal", traitm="bm", matdir=matdir)
			
			## OU
			sim_traits_nets_par(all_t1_bal_ou, method = "b", type="bal", traitm="ou", matdir=matdir)
			sim_traits_nets_par(all_t1_unbal_ou, method = "b", type="unbal", traitm="ou", matdir=matdir)
			
			## EB
			sim_traits_nets_par(all_t1_bal_eb, method = "b", type="bal", traitm="eb", matdir=matdir)
			sim_traits_nets_par(all_t1_unbal_eb, method = "b", type="unbal", traitm="eb", matdir=matdir)
		}
	
# 	################## Calculate network metrics on matrices
# 	message("...calculating network structures...")
# 	df_rand <- getnetmets(mats_rand_bal, mats_rand_unbal, netmets=netmets) # random networks
# 	
# 	## BM
# 	df_traits1_bm_comp <- getnetmets(mats_traits1_bal_bm_comp, mats_traits1_unbal_bm_comp, netmets=netmets) # complementarity traits
# 	df_traits1_bm_barr <- getnetmets(mats_traits1_bal_bm_barr, mats_traits1_unbal_bm_barr, netmets=netmets) # barrier traits
# 	
# 	## OU
# 	df_traits1_ou_comp <- getnetmets(mats_traits1_bal_ou_comp, mats_traits1_unbal_ou_comp, netmets=netmets) # complementarity traits
# 	df_traits1_ou_barr <- getnetmets(mats_traits1_bal_ou_barr, mats_traits1_unbal_ou_barr, netmets=netmets) # barrier traits
# 	
# 	## EB
# 	df_traits1_eb_comp <- getnetmets(mats_traits1_bal_eb_comp, mats_traits1_unbal_eb_comp, netmets=netmets) # complementarity traits
# 	df_traits1_eb_barr <- getnetmets(mats_traits1_bal_eb_barr, mats_traits1_unbal_eb_barr, netmets=netmets) # barrier traits
# 	
	################## Writing matrices to files if dumpmatrices=TRUE
	
# 	if(dumpmatrices==TRUE){
# 		message("...writing matrices to files...")
# 		matsall <- list(mats_rand_bal,mats_rand_unbal,mats_traits1_bal_bm_comp,mats_traits1_unbal_bm_comp,
# 										mats_traits1_bal_ou_comp,mats_traits1_unbal_ou_comp,mats_traits1_bal_eb_comp,
# 										mats_traits1_unbal_eb_comp,mats_traits1_bal_bm_barr,mats_traits1_unbal_bm_barr,
# 										mats_traits1_bal_ou_barr,mats_traits1_unbal_ou_barr,mats_traits1_bal_eb_barr,
# 										mats_traits1_unbal_eb_barr)
# 		matsnames <- list("mats_rand_bal","mats_rand_unbal","mats_traits1_bal_bm_comp",
# 											"mats_traits1_unbal_bm_comp","mats_traits1_bal_ou_comp","mats_traits1_unbal_ou_comp",
# 											"mats_traits1_bal_eb_comp","mats_traits1_unbal_eb_comp","mats_traits1_bal_bm_barr",
# 											"mats_traits1_unbal_bm_barr","mats_traits1_bal_ou_barr","mats_traits1_unbal_ou_barr",
# 											"mats_traits1_bal_eb_barr","mats_traits1_unbal_eb_barr")
# 		names(matsall) <- matsnames
# 		# 		dir.create(file.path(matdir), showWarnings=FALSE)
# 		for(i in 1:length(matsnames)) {
# 			for(j in 1:length(matsall[[i]])) {
# 				write.table(matsall[[i]][[j]], file=paste(matdir, "/", "sp", sum(tips_p,tips_a), "_", matsnames[[i]], "_", j, ".web", sep=""), row.names=F, col.names=F)
# 				# 				write.table(matsall[[i]][[j]], file=paste(matsnames[[i]], "_", j, ".web", sep=""), row.names=F, col.names=F)
# 			}
# 		}
# 	} else
# 	{ NULL }	
# 	
	################## Combining data to a data.frame
# 	message("...combining data to a data.frame...")
# 	
# 	alldat <- rbind(df_rand, 
# 									df_traits1_bm_comp, df_traits1_bm_barr,
# 									df_traits1_ou_comp, df_traits1_ou_barr,
# 									df_traits1_eb_comp, df_traits1_eb_barr)
# 	alldat$traitevol <- c( rep("Random",nrow(df_rand)), 
# 												 rep("BM",nrow(df_traits1_bm_comp)),
# 												 rep("BM",nrow(df_traits1_bm_barr)),
# 												 rep("OU",nrow(df_traits1_ou_comp)),
# 												 rep("OU",nrow(df_traits1_ou_barr)),
# 												 rep("EB",nrow(df_traits1_eb_comp)),
# 												 rep("EB",nrow(df_traits1_eb_barr))
# 	)
# 	alldat$netmethod <- c( rep("Random",nrow(df_rand)), 
# 												 rep("Complementarity",nrow(df_traits1_bm_comp)),
# 												 rep("Barrier",nrow(df_traits1_bm_barr)),
# 												 rep("Complementarity",nrow(df_traits1_ou_comp)),
# 												 rep("Barrier",nrow(df_traits1_ou_barr)),
# 												 rep("Complementarity",nrow(df_traits1_eb_comp)),
# 												 rep("Barrier",nrow(df_traits1_eb_barr))
# 	)
# 	alldat$numsp_p <- rep(tips_p, nrow(alldat))
# 	alldat$numsp_a <- rep(tips_a, nrow(alldat))
# 	alldat$numsp_all <- rep(sum(tips_a,tips_p), nrow(alldat))
	write.csv(traits_df, paste(matdir, "/sp", sum(tips_p,tips_a), substring(modeltorun, 1, 4), ".csv", sep=""), row.names=F)
	message("...done.")
# 	return( traits_df )
}