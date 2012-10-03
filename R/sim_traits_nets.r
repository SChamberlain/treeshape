# Simulate networks, with interactions propoprtional to trait matching
sim_traits_nets <- function(listoftraitvecs, method = c("ratio","complementarity","barrier"), value) {
  mats <- list()
  method <- match.arg(method, c("ratio","complementarity","barrier"))
  message(paste("Using the ", method, " method"))
  for(i in 1:length(listoftraitvecs[[1]])) {
    # where the interaction occurs or not
    ## Ratio - e.g., body size ratio, for gape limitation
    if(method == "ratio"){
      mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]], 
                  function(x,y) as.numeric(exp(x-y) < value)) 
    } else
      if(method == "complementarity"){
        mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]], 
                    function(x,y) as.numeric(abs(x-y) < value))
      }  else
        if(method == "barrier"){
          mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]], 
                      function(x,y) as.numeric(x > y))
        } else 
          stop("must be one of ratio, complementarity or barrier")
    dimnames(mm)[[1]] <- names(listoftraitvecs[[1]][[i]])
    dimnames(mm)[[2]] <- names(listoftraitvecs[[2]][[i]])
    if(sum(mm) == 0) { mm <- NULL } else 
      if( sum(mm) == nrow(mm) * ncol(mm) ) {mm <- NULL } else
      { mm <- mm }
    mats[[i]] <- mm
  }
  mats[!sapply(mats, is.null)]
}