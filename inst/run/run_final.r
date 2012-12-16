# library(devtools); install_github('treeshape', 'schamberlain')
library(doMC); library(plyr); library(ape); library(apTreeshape); library(bipartite); 
library(ggplot2); library(reshape2); library(picante); library(geiger); library(parallel)
library(treeshape)

numtrees <- 1000
# numspvec <- c(15, 20, 30, 40, 50, 60, 70, 80, 100, 110, 130, 150, 175, 200)

###### run it
# colless one at a time
out_colless_15_b <- simbaltrees_topy(15, metric="colless", numtrees=numtrees, cutlow=-0.7, cuthigh=0.7, a=10, bounds=c(0,100), alpha=10, sigma=3,
																		 alpha_eb=-1.1, sigma_eb=3, cval=0.5, asymm=2, matdir="~/newfiles2", modeltorun="barrier")
out_colless_15_c <- simbaltrees_topy(15, metric="colless", numtrees=numtrees, cutlow=-0.7, cuthigh=0.7, a=10, bounds=c(0,100), alpha=10, sigma=3,
																		 alpha_eb=-1.1, sigma_eb=3, cval=0.5, asymm=2, matdir="~/newout_dec/colless_mats/15", modeltorun="complementarity")
out_colless_20 <- simbaltrees_topy(20, metric="colless", numtrees=numtrees, cutlow=-0.7, cuthigh=0.7, a=10, bounds=c(0,100), alpha=10, sigma=3,
																	 alpha_eb=-1.1, sigma_eb=3, cval=0.5, asymm=2, matdir="~/newout_dec/colless_mats/20")
out_colless_30 <- simbaltrees_topy(30, metric="colless", numtrees=numtrees, cutlow=-0.7, cuthigh=0.7, a=10, bounds=c(0,100), alpha=10, sigma=3,
																	 alpha_eb=-1.1, sigma_eb=3, cval=0.5, asymm=2, matdir="~/newout_dec/colless_mats/30")
out_colless_40 <- simbaltrees_topy(40, metric="colless", numtrees=numtrees, cutlow=-0.7, cuthigh=0.7, a=10, bounds=c(0,100), alpha=10, sigma=3,
																	 alpha_eb=-1.1, sigma_eb=3, cval=0.5, asymm=2, matdir="~/newout_dec/colless_mats/40")
out_colless_50 <- simbaltrees_topy(50, metric="colless", numtrees=numtrees, cutlow=-0.7, cuthigh=0.7, a=10, bounds=c(0,100), alpha=10, sigma=3,
																	 alpha_eb=-1.1, sigma_eb=3, cval=0.5, asymm=2, matdir="~/newout_dec/colless_mats/50")


install.packages(c("doMC", "plyr", "ape", "apTreeshape", "bipartite", "ggplot2", "reshape2", 
									 "picante", "geiger"))
install_github('treeshape', 'schamberlain')
library(doMC); library(plyr); library(ape); library(apTreeshape); library(bipartite); 
library(ggplot2); library(reshape2); library(picante); library(geiger); 
library(treeshape)