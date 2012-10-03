################## Calculate network metrics on matrices
getnetmets <- function(balanced, unbalanced) {
  netmets_bal <- ldply(balanced, 
                  function(x) networklevel(x, 
                  index = c("connectance", "links per species", "nestedness", "web asymmetry")))
  netmets_unbal <- ldply(unbalanced, 
                  function(x) networklevel(x, 
                  index = c("connectance", "links per species", "nestedness", "web asymmetry")))
  data.frame( 
    type = c( rep("bal", length(balanced)), rep("unbal", length(unbalanced))), 
    rbind(netmets_bal, netmets_unbal) )
}