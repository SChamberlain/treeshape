#' Random traits - randomize traits across tips - shuffles names, but doesn't 
#' 		move position in tree of course
#' 
#' This is used in randomizing traits across the tips of a phylogeny so as to
#' 		randomize any phylogenetic signal.
#' 
#' @param w A named vector of traits. Each element in vector must be named.
#' @examples \dontrun{
#' tree <- pbtree(n=10)
#' traitvec <- fastBM(tree)
#' randtips(traitvec)
#' }
#' @export
randtips <- function(w){
  names_ <- names(w)
  names_shuff <- sample(names_, size=length(names_), replace=F)
  names(w) <- names_shuff
  w
}