#' Measure aspects of traits on trees.
#' 
#' @import picante
#' @param x A named trait vector, the names should match the tip names in the tree.
#' @param y A phylogeny, has to be a phylo object.
#' @examples \dontrun{
#' trees <- rmtree(N=10, n=10)
#' traitvecs <- lapply(trees, fastBM)
#' traitsig(traitvecs, trees)
#' }
#' @export
traitsig <- function(x, y) {
  out <- data.frame(psig = NA, traitvar = NA)
  for(i in 1:length(x)) {
    out[i,] <- c(Kcalc(x[[i]], y[[i]]), var(x[[i]]))
  }
  out
}