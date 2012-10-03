################## Measure aspects of traits on trees
traitsig <- function(x, y) {
  out <- data.frame(psig = NA, traitvar = NA)
  for(i in 1:length(x)) {
    out[i,] <- c(Kcalc(x[[i]], y[[i]]), var(x[[i]]))
  }
  out
}