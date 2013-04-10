#' Calculate network metrics on matrices.
#' 
#' This is used in calculating network metrics on pairs of matrices made from
#' 		balanced and unbalanced trees.
#' 
#' @import plyr bipartite stringr
#' @param filepath Path to an edgelist
#' @param netmets Network structure metrics to calculate - only use those that 
#' 		calculate single values for each matrix.
#' @examples \dontrun{
#' getnetmets2(filepath = "/Users/scottmac2/testtest/sp30_bal_eb_twom_4.web")
#' 
#' filepaths <- c("/Users/scottmac2/testtest/sp30_bal_eb_twom_4.web", 
#' 								"/Users/scottmac2/testtest/sp30_bal_bm_twom_3.web",
#' 								"/Users/scottmac2/testtest/sp30_bal_bm_twom_2.web",
#' 								"/Users/scottmac2/testtest/sp30_bal_bm_twom_1.web")
#' l_ply(filepaths, getnetmets2)
#' }
#' @export
getnetmets2 <- function(filepath, netmets=c("connectance", "NODF2"), directory)
{
	webname <- str_split(filepath, "/")[[1]][length(str_split(filepath, "/")[[1]])]
	webname_split <- str_split(webname, "_")[[1]]
	number <- str_split(webname_split[length(webname_split)], "\\.")[[1]][[1]]
	fileinfo <- c(webname_split[-length(webname_split)], number)
	
	mat <- edgelist2matrix(filepath)
	if( any(grepl("NODF2", netmets, ignore.case=T)) == TRUE )
		nodf2_ <- nested(mat, method="NODF2")
	
	netmets_other <- netmets[!grepl("NODF2", netmets, ignore.case=T)]
	netmets_ <- networklevel(mat, index = netmets_other)
  
	temp <- data.frame(t(data.frame(c(fileinfo, nodf2_[[1]], netmets_[[1]]))))
	row.names(temp) <- NULL
	
	write.table(temp, file=paste0(directory, "/netmets.csv"), sep = ',', 
							append = TRUE, row.names = FALSE, col.names = FALSE)
}