#' Simulate networks, with interactions propoprtional to trait matching
#' 
#' @import parallel
#' @param listoftraitvecs 
#' @param method The model to be used to construct interaction matrices. One of
#' 		"ratio","complementarity","barrier".
#' @param value Value at which to determine if species interact or not - 
#' 		this value depends on the model you are using.
#' @return A data.frame of network structure metrics for balanced and unbalanced 
#' 		trees.
#' @examples \dontrun{
#' trees <- rmtree(N=1000, n=60)
#' trees2 <- rmtree(N=1000, n=120)
#' registerDoMC(4)
#' traitvecs <- llply(trees, fastBM, bounds=c(0,Inf), .parallel=TRUE)
#' traitvecs2 <- llply(trees2, fastBM, bounds=c(0,Inf), .parallel=TRUE)
#' alltraits <- list(traitvecs, traitvecs2)
#' out<-sim_traits_nets_par(listoftraitvecs=alltraits, method="c", value=0.5, cores=4)
#' length(out)
#' }
#' @export
sim_traits_nets_par <- function(listoftraitvecs, cores = 2, 
														 method = c("ratio","complementarity","barrier"), value) 
{
# 	mats <- list()
	method <- match.arg(method, c("ratio","complementarity","barrier"))
	
# 	system.time( mapply(funcy, listoftraitvecs[[1]], listoftraitvecs[[2]], SIMPLIFY=F) )
	
	funcy <- function(list1, list2) {
		if(method == "complementarity"){
			
			mmm <- NULL
			while(is.null(mmm)){
				
				mm <- outer(list1, list2, 
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
					
					mm <- outer(list1, list2, 
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
		return( doit(mmm) )
	} # end of function	

	mats <- mcmapply(funcy, listoftraitvecs[[1]], listoftraitvecs[[2]], SIMPLIFY=F, mc.cores = cores) 
	mats[!sapply(mats, is.null)]
}