#' Measure aspects of traits on trees.
#' 
#' @import ape
#' @param listoftrees 
#' @examples \dontrun{
#' trees <- list(rmtree(N=50, n=10), rmtree(N=50, n=20))
#' sim_rand_nets(trees)
#' }
#' @export
sim_rand_nets <- function(listoftrees) {
  mats <- list()
  for(i in 1:length(listoftrees[[1]])) {
    m <- Ntip(listoftrees[[1]][[i]]) # number of plant species
    n <- Ntip(listoftrees[[2]][[i]]) # number of animal species
    # make random matrix and put matrix into list
    mm <- matrix(rbinom(m * n, 1, .5), ncol = m, nrow = n)  
    dimnames(mm)[[1]] <- as.character(listoftrees[[2]][[i]]$tip.label)
    dimnames(mm)[[2]] <- as.character(listoftrees[[1]][[i]]$tip.label)
    mats[[i]] <- mm
  }
  mats
}