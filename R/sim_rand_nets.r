#' Measure aspects of traits on trees.
#' 
#' @import ape
#' @param listoftrees 
#' @examples \dontrun{
#' trees <- rmtree(N=10, n=10)
#' traitvecs <- lapply(trees, fastBM)
#' traitsig(traitvecs, trees)
#' }
#' @export
sim_rand_nets <- function(listoftrees) {
  mats <- list()
  for(i in 1:length(listoftrees[[1]])) {
    m <- Ntip(listoftrees[[1]][[i]]) # number of plant species
    n <- Ntip(listoftrees[[2]][[i]]) # number of animal species
    # make random matrix and put matrix into list
    mm <- matrix(rbinom(m * n, 1, .5), ncol = m, nrow = n)  
    dimnames(mm)[[1]] <- as.list(listoftrees[[1]][[i]]$tip.label)
    dimnames(mm)[[2]] <- as.list(listoftrees[[2]][[i]]$tip.label)
    mats[[i]] <- mm
  }
  mats
}