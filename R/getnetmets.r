#' Calculate network metrics on matrices.
#' 
#' This is used in calculating network metrics on pairs of matrices made from
#' 		balanced and unbalanced trees.
#' 
#' @import plyr bipartite
#' @param balanced A list of matrices made from balanced trees.
#' @param unbalanced A list of matrices made from unbalanced trees.
#' @param netmets Network structure metrics to calculate - only use those that 
#' 		calculate single values for each matrix.
#' @examples \dontrun{
#' # Let's pretend these are balanced and unbalanced matrices
#' m <- 10
#' n <- 5
#' netmets <- c("connectance", "links per species", "nestedness", "web asymmetry")
#' mats_rand_bal <- replicate(20, matrix(rbinom(m * n, 1, .5), ncol = m, nrow = n), FALSE)
#' mats_rand_unbal <- replicate(20, matrix(rbinom(m * n, 1, .5), ncol = m, nrow = n), FALSE)
#' getnetmets(balanced = mats_rand_bal, mats_rand_unbal, netmets)
#' }
#' @export
getnetmets <- function(balanced, unbalanced, netmets) 
{
  netmets_bal <- ldply(balanced, function(x) networklevel(x, index = netmets))
  netmets_unbal <- ldply(unbalanced, function(x) networklevel(x, index = netmets))
  data.frame( 
    type = c( rep("bal", length(balanced)), rep("unbal", length(unbalanced))), 
    rbind(netmets_bal, netmets_unbal) )
}