#' Simulate networks, with interactions propoprtional to trait matching
#' 
#' @import parallel bipartite
#' @param listoftraitvecs Nested list of trait vectors.
#' @param type Balanced ("bal") or unbalanced ("unbal"). 
#' @param traitm Evolutionary model, used to create file name, one of "bm", "ou", or "eb".
#' @param matdir directory to put matrices in
#' @param method The model to be used to construct interaction matrices. One of
#' 		"complementarity", "barrier", or "twomethods" (for combined complementarity/barrier).
#' @param includetraitvar Include trait variation in the simulation or not. 
#' 		Defaults to FALSE.
#' @param value Value at which to determine if species interact or not - 
#' 		this value depends on the model you are using.
#' @param output One of matrix or edgelist.
#' @return A data.frame of network structure metrics for balanced and unbalanced 
#' 		trees (if output="matrix"), or an edgelist (if output="edgelist").
#' @examples \dontrun{
#' trees <- rmtree(N=100, n=15)
#' trees2 <- rmtree(N=100, n=30)
#' library(doMC)
#' registerDoMC(4)
#' traitvecs <- llply(trees, fastBM, bounds=c(0,Inf), .parallel=TRUE)
#' traitvecs2 <- llply(trees2, fastBM, bounds=c(0,Inf), .parallel=TRUE)
#' alltraits <- list(traitvecs, traitvecs2)
#' sim_traits_nets_par(listoftraitvecs=alltraits, method="c", value=0.5, type="bal", traitm="bm", matdir="~/testtest", includetraitvar = TRUE)
#' sim_traits_nets_par(listoftraitvecs=alltraits, method="b", type="bal", traitm="bm", matdir="~/newfiles2")
#' length(out)
#' }
#' @export
sim_traits_nets_par <- function(listoftraitvecs, type = NULL, traitm = NULL, matdir = "~",
		method = c("complementarity","barrier","twomethods"), includetraitvar = FALSE,
		value = NULL, output = "matrix") 
{
	doit <- function(x) {
		ui <- function(x) if( sum(x)==0 ){replace(x, sample(grep(0, x), 1), 1)} else{x}
		if(is.null(x)){NULL} else
		{	
			xx <- apply(x, 2, ui)
			t(apply(xx, 1, ui))
		}
	}
	
	method <- match.arg(method, c("complementarity","barrier","twomethods"))
	
	for(i in 1:length(listoftraitvecs[[1]])) {
		if(method == "complementarity"){
			
			if(includetraitvar){
				makevar <- function(listofsp){
					# make trait variation values, with min of 0 and max of 0.25
					# matches the "narrow-range" complementarity model of Santamaria & Rodriguez PLoS paper
					temp <- runif(n=length(listofsp),min=0,max=0.25) 
					names(temp) <- names(listofsp)
					return(temp)
				}
				one <- listoftraitvecs[[1]][[i]]
				two <- listoftraitvecs[[2]][[i]]
				one_var <- makevar(one)
				two_var <- makevar(two)
				aaa <- list(mean=one, var=one_var)
				bbb <- list(mean=two, var=two_var)
				todo <- list(aaa,bbb)
				withvar <- function(i,j,data){ as.numeric(abs(data[[1]][[1]][i]-data[[2]][[1]][j]) < value*(data[[1]][[2]][i]+data[[2]][[2]][j]) ) }
				withvar_ <- Vectorize(withvar, vectorize.args=list("i","j"))
				mm <- outer(names(one),names(two), withvar_, data=todo)
			} else
			{
				one <- listoftraitvecs[[1]][[i]]
				two <- listoftraitvecs[[2]][[i]]
				mm <- outer(one, two, function(x,y) as.numeric(abs(x-y) < value))
			}
			
			# Remove any matrices that have all zeros or all ones
			if(sum(mm) == 0) { mmm <- NULL } else 
				if( sum(mm) == nrow(mm) * ncol(mm) ) {mmm <- NULL } else
				{ mmm <- mm }
			if(is.null(mmm)){NULL} else 
			{
				mmm <- doit(mmm)
				
				if(output=="matrix"){
					write.table(mmm, file=paste(matdir, "/", "sp", sum(length(one),length(two)),
																			"_",type,"_",traitm,"_",substring(method, 1, 4),"_",i,".web", sep=""), row.names=F, col.names=F)
				} else
				{
					edgelist <- web2edges(mmm, out.files="", weight.column=FALSE, return=TRUE)
					write.table(edgelist, row.names=FALSE, col.names=FALSE, file=paste(matdir, "/", "sp", sum(length(one),length(two)),
																														"_",type,"_",traitm,"_",substring(method, 1, 4),"_",i,".web", sep=""))
				}
			}
		}  else
			if(method == "barrier"){			
				one <- listoftraitvecs[[1]][[i]]
				two <- listoftraitvecs[[2]][[i]]
				mm <- outer(one, two, function(x,y) as.numeric(x < y))
				# Remove any matrices that have all zeros or all ones
				if(sum(mm) == 0) { mmm <- NULL } else 
					if( sum(mm) == nrow(mm) * ncol(mm) ) {mmm <- NULL } else
					{ mmm <- mm }
				if(is.null(mmm)){NULL} else 
				{
					mmm <- doit(mmm)
					
					if(output=="matrix"){
					write.table(mmm, file=paste(matdir, "/", "sp", sum(length(one),length(two)),
																			"_",type,"_",traitm,"_",substring(method, 1, 4),"_",i,".web", sep=""), row.names=F, col.names=F)
					} else
					{
						edgelist <- web2edges(mmm, out.files="", weight.column=FALSE, return=TRUE)
						write.table(edgelist, row.names=FALSE, col.names=FALSE, file=paste(matdir, "/", "sp", sum(length(one),length(two)),
																															"_",type,"_",traitm,"_",substring(method, 1, 4),"_",i,".web", sep=""))
					}
				}
			} else
				if(method == "twomethods"){
					
					if(includetraitvar){
						makevar <- function(listofsp){
							# make trait variation values, with min of 0 and max of 0.25
							# matches the "narrow-range" complementarity model of Santamaria & Rodriguez PLoS paper
							temp <- runif(n=length(listofsp),min=0,max=0.25) 
							names(temp) <- names(listofsp)
							return(temp)
						}
						one <- listoftraitvecs[[1]][[i]]
						two <- listoftraitvecs[[2]][[i]]
						one_var <- makevar(one)
						two_var <- makevar(two)
						aaa <- list(mean=one, var=one_var)
						bbb <- list(mean=two, var=two_var)
						todo <- list(aaa,bbb)
						cb_withvar <- function(i,j,data){ as.numeric((abs(data[[1]][[1]][i]-data[[2]][[1]][j]) < value*(data[[1]][[2]][i]+data[[2]][[2]][j])) & data[[1]][[1]][i] < data[[2]][[1]][j]) }
						cb_withvar_ <- Vectorize(cb_withvar, vectorize.args=list("i","j"))
						mm <- outer(names(one),names(two), cb_withvar_, data=todo)  # complementarity-barrier combined model
					} else
					{
						one <- listoftraitvecs[[1]][[i]]
						two <- listoftraitvecs[[2]][[i]]
						mm <- outer(one, two, function(x,y) as.numeric((abs(x-y) < value) & x < y))
					}
					
					# Remove any matrices that have all zeros or all ones
					if(sum(mm) == 0) { mmm <- NULL } else 
						if( sum(mm) == nrow(mm) * ncol(mm) ) {mmm <- NULL } else
						{ mmm <- mm }
					if(is.null(mmm)){NULL} else 
					{
						mmm <- doit(mmm)
						
						if(output=="matrix"){
							write.table(mmm, file=paste(matdir, "/", "sp", sum(length(one),length(two)),
																					"_",type,"_",traitm,"_",substring(method, 1, 4),"_",i,".web", sep=""), row.names=F, col.names=F)
						} else
						{
							edgelist <- web2edges(mmm, out.files="", weight.column=FALSE, return=TRUE)
							write.table(edgelist, row.names=FALSE, col.names=FALSE, file=paste(matdir, "/", "sp", sum(length(one),length(two)),
																																								 "_",type,"_",traitm,"_",substring(method, 1, 4),"_",i,".web", sep=""))
						}
					}
				} else
				stop("must be one of complementarity, barrier, twomethods, or their abbreviations")
	} # end of loop
}