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
#' trees <- rmtree(N=20, n=5)
#' trees2 <- rmtree(N=20, n=10)
#' traitvecs <- lapply(trees, fastBM, bounds=c(0,Inf))
#' traitvecs2 <- lapply(trees2, fastBM, bounds=c(0,Inf))
#' alltraits <- list(traitvecs, traitvecs2)
#' sim_traits_nets2(listoftraitvecs=alltraits, method="c", value=1)
#' }
#' @export
sim_traits_nets2 <- function(listoftraitvecs, 
	method = c("ratio","complementarity","barrier"), value) 
{
	mats <- list()
	method <- match.arg(method, c("ratio","complementarity","barrier"))
	
	
	for(i in 1:length(listoftraitvecs[[1]])) { # start loop
		# where the interaction occurs or not
		## Ratio - e.g., body size ratio, for gape limitation
		if(method == "ratio"){
			mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]],
									function(x,y) as.numeric(exp(x-y) < value)) 
		} else
			if(method == "complementarity"){
				
				mmm <- NULL
				while(is.null(mmm)){
					
					mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]], 
											function(x,y) as.numeric(abs(x-y) < value))
					# Remove any matrices that have all zeros or all ones
					if(sum(mm) == 0) { mmm <- NULL } else 
						if( sum(mm) == nrow(mm) * ncol(mm) ) {mmm <- NULL } else
							{ mmm <- mm }
				}
				
			}  else
				if(method == "barrier"){
					
					mmm <- NULL
					while(is.null(mmm)){
						
						mm <- outer(listoftraitvecs[[1]][[i]], listoftraitvecs[[2]][[i]], 
												function(x,y) as.numeric(x > y))
						# Remove any matrices that have all zeros or all ones
						if(sum(mm) == 0) { mmm <- NULL } else 
							if( sum(mm) == nrow(mm) * ncol(mm) ) {mmm <- NULL } else
							{ mmm <- mm }
					}
					
				} else 
					stop("must be one of ratio, complementarity or barrier, or their abbreviations")
		
		# Add a random 1 to a species that has all zeros
		doit <- function(x) {
			ui <- function(x) if( sum(x)==0 ){replace(x, sample(grep(0, x), 1), 1)} else{x}
			if(is.null(x)){NULL} else
			{	
				xx <- apply(x, 2, ui)
				t(apply(xx, 1, ui))
			}
		}
		mm <- doit(mmm)
		mats[[i]] <- mm
	} # end loop
	
	mats[!sapply(mats, is.null)]
}