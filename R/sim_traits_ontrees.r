#' Simulate traits on trees with various evolutionary models.
#' 
#' Use any of Brownian Motion, Orntsein-Uhlenbeck, or Early-Burst. The function
#' 		uses functions already created for evolving traits on trees from the 
#' 		packages: ape, phytools, and geiger. 
#' 
#' @import ape geiger
#' @param trees A list of phylogenetic trees.  
#' @param model The model to be used to evolve traiats on a phylogenetic tree. 
#' 		One of "bm","ou", or "eb" for Brownian Motion, Orntsein-Uhlenbeck, and 
#' 		Early-Burst, respectively. "bm" uses the fastBM function from phytools. 
#' 		"ou" uses the rTraitCont function from ape. "eb" transforms a tree using 
#' 		exponentialchangeTree function from geiger, then evolves traits with 
#' 		rTraitCont function from ape.
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
#' @return A list of named vectors of traits matching the tips of the input trees.
#' @examples \dontrun{
#' sim_traits_ontrees(trees = trees_colless_plants_bal, model = "bm")
#' sim_traits_ontrees(trees = trees_colless_plants_bal, model = "ou", alpha=1.5, sigma=2)
#' sim_traits_ontrees(trees = trees_colless_plants_bal, model = "eb", alpha_eb=-0.9, sigma_eb=3)
#' }
#' @export
sim_traits_ontrees <- function(trees, model = c("bm","ou","eb"), a = 10,
	bounds = c(0,100), alpha=1, sigma=1, theta=1, alpha_eb=-0.8, sigma_eb=3) 
{
	model <- match.arg(model, c("bm","ou","eb")) 
	if(model=="bm"){
		lapply(trees, fastBM, a = 10, bounds=c(0,100))
	} else
		if(model=="ou"){
			temp <- lapply(trees, rTraitCont, model = "OU", sigma = sigma, alpha=alpha, theta=theta)
			lapply(temp, function(x) x + abs(min(x)))
		} else
			if(model=="eb"){
				temp <- lapply(trees, function(x) rTraitCont(exponentialchangeTree(x, a=alpha_eb), sigma = sigma_eb))
				lapply(temp, function(x) x + abs(min(x)))
			} else
				 {stop("need to specify model as either bm, ou, or eb")}
}