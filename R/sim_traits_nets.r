#' Simulate networks, with interactions propoprtional to trait matching
#' 		
#' @param listoftraitvecs 
#' @param method The model to be used to construct interaction matrices. One of
#' 		"ratio","complementarity","barrier".
#' @param value Value at which to determine if species interact or not - 
#' 		this value depends on the model you are using.
#' @return A data.frame of network structure metrics for balanced and unbalanced 
#' 		trees.
#' @examples \dontrun{
#' trees <- rmtree(N=10, n=10)
#' trees2 <- rmtree(N=10, n=10)
#' traitvecs <- lapply(trees, fastBM)
#' traitvecs2 <- lapply(trees, fastBM)
#' alltraits <- list(traitvecs, traitvecs2)
#' sim_traits_nets(alltraits, method="r", value=0.5) # where r = ratio, you can abbreviate
#' }
#' @export
sim_traits_nets <- function(listoftraitvecs, 
			method = c("ratio","complementarity","barrier"), value) 
{
  mats <- list()
  method <- match.arg(method, c("ratio","complementarity","barrier"))
  message(paste("Using the ", method, " method"))
  for(i in 1:length(listoftraitvecs[[1]])) {
    if(method == "ratio"){
      mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]], 
      						function(x,y) as.numeric(y/x < value)) 
    } else
      if(method == "complementarity"){
        mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]], 
                    function(x,y) as.numeric(abs(x-y) < value))
      }  else
        if(method == "barrier"){
          mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]], 
                      function(x,y) as.numeric(y > x))
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