#' Plotting output from simbaltrees
#' 
#' @param output Output from simbaltrees.
#' @param netval A network metric used in your call to simbaltrees.
#' @param justdf Return just the data.frame or not (default to FALSE).
#' @examples \dontrun{
#' netmets <- c("connectance", "nestedness","nodf2")
#' outout <- suppressMessages(llply(c(10, 15), function(x) 
#' 		simbaltrees(tips_p=10, metric="colless", numtrees=5, cutlow=-0.5, cuthigh=0.5, 
#' 		a=10, bounds=c(0,100), alpha=1, sigma=1, alpha_eb=-0.8, sigma_eb=3, cval=0.8, 
#' 		asymm=2, dumpmatrices=TRUE, matdir="~/newfiles2", netmets=netmets)
#' simbal_plot(outout, "connectance")
#' }
simbal_plot <- function(output, netval, justdf=FALSE)
{
	outdf <- ldply(output, function(x) x[[2]])
	outdf_melt <- melt(outdf, id.vars=c(1,6:8))  
	outdf_melt_ <- ddply(outdf_melt, .(type, traitevol, netmethod, numsp, variable), summarise, 
											 mean = mean(value, na.rm=T),
											 se = sd(value, na.rm=T)/sqrt(na.omit(length(value))),
											 n_trees = length(na.omit(value))
	)
	outdf_melt_$conf_int <- outdf_melt_$se*1.96
	dfplot <- outdf_melt_
	limits <- aes(ymax = mean + conf_int, ymin = mean - conf_int)
	dodge <- position_dodge(w=0.3)
	dfplot$traitevol <- factor(dfplot$traitevol, levels=c("Random","BM","EB"))
	dfplot <- dfplot[!dfplot$netmethod %in% "Random",]
	
	if(justdf == TRUE){ dfplot } else { 
		
		ggplot(dfplot[dfplot$variable==netval,], aes(numsp, mean, shape = type, colour = type)) +
			geom_point(size=4, position=dodge) +
			geom_errorbar(limits, width=0.2, position=dodge) +
			theme_bw() +
			facet_grid(traitevol ~ netmethod, scales="free")
		
	}
}
# 
# simbal_plotb(rval01_cval11, "connectance")
# simbal_plotb(rva05_cval15, "connectance")
# simbal_plotb(rval1_cval2, "connectance")
# simbal_plotb(rval15_cval25, "connectance")
# simbal_plotb(rval2_cval3, "connectance")
# simbal_plotb(rval2_cval3, "connectance")
# simbal_plotb(rval25_cval5, "connectance")
# simbal_plotb(rval3_cval8, "connectance")
# simbal_plotb(rval5_cval8, "connectance")
# 
# simbal_plotb(rval01_cval11, "nestedness")
# simbal_plotb(rva05_cval15, "nestedness")
# simbal_plotb(rval1_cval2, "nestedness")
# simbal_plotb(rval15_cval25, "nestedness")
# simbal_plotb(rval2_cval3, "nestedness")
# simbal_plotb(rval2_cval3, "nestedness")
# simbal_plotb(rval25_cval5, "nestedness")
# simbal_plotb(rval3_cval8, "nestedness")
# simbal_plotb(rval5_cval8, "nestedness")
# 
# simbal_plotb(rval01_cval11, "web.asymmetry")
# simbal_plotb(rva05_cval15, "web.asymmetry")
# simbal_plotb(rval1_cval2, "web.asymmetry")
# simbal_plotb(rval15_cval25, "web.asymmetry")
# simbal_plotb(rval2_cval3, "web.asymmetry")
# simbal_plotb(rval2_cval3, "web.asymmetry")
# simbal_plotb(rval25_cval5, "web.asymmetry")
# simbal_plotb(rval3_cval8, "web.asymmetry")
# simbal_plotb(rval5_cval8, "web.asymmetry")
# 
# 
# setwd("newoutput")
# write.csv(simbal_plotb(rval01_cval11, justdf=TRUE), "rval01_cval11.csv")
# write.csv(simbal_plotb(rva05_cval15, justdf=TRUE), "rva05_cval15.csv")
# write.csv(simbal_plotb(rval1_cval2, justdf=TRUE), "rval1_cval2.csv")
# write.csv(simbal_plotb(rval15_cval25, justdf=TRUE), "rval15_cval25.csv")
# write.csv(simbal_plotb(rval2_cval3, justdf=TRUE), "rval2_cval3.csv")
# write.csv(simbal_plotb(rval25_cval5, justdf=TRUE), "rval25_cval5.csv")
# write.csv(simbal_plotb(rval3_cval8, justdf=TRUE), "rval3_cval8.csv")