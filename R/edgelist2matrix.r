#' Make a matrix from an edgelist on file.
#' 
#' @param path Path to a edgelist file.
#' @examples
#' edgelist2matrix("/Users/scottmac2/testtest/sp30_bal_ou_comp_1.web") 
#' @export
edgelist2matrix <- function(path){
	x <- read.table(file=path)
	emptymat <- matrix(0,nrow=length(unique(x[,1])),ncol=length(unique(x[,2])),dimnames=list(unique(x[,1]),unique(x[,2])))
	for(i in 1:nrow(x)){
		emptymat[grep(x[i,1], dimnames(emptymat)[[1]], value=TRUE), grep(x[i,2], dimnames(emptymat)[[2]], value=TRUE)] <- 1
	}
	return(emptymat)
}