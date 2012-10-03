## Random traits - randomize traits across tips - shuffles names, but doesn't move position in tree of course
randtips <- function(w){
  names_ <- names(w)
  names_shuff <- sample(names_, size=length(names_), replace=F)
  names(w) <- names_shuff
  w
}