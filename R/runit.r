# install.packages("doMC")
library(doMC)
library(plyr); library(ape); library(apTreeshape); library(bipartite); library(ggplot2); library(reshape2); library(picante)
source("fastBM.r"); source("nodeHeights.r"); source("pbtree.r"); source("sim_rand_nets.r"); source("sim_traits_nets.r"); source("getnetmets.r"); source("randtips.r"); source("traitsig.r")
registerDoMC(cores=4)
# out <- suppressMessages(llply(c(10,15), function(x) simbaltrees(x, metric="colless", numtrees=10, cutlow=-0.5, cuthigh=0.5), .parallel=TRUE))
out <- suppressMessages(llply(c(20, 40, 60, 80, 100), function(x) simbaltrees(x, metric="colless", numtrees=50, cutlow=-0.5, cuthigh=0.5), .parallel=TRUE))
out_larger <- suppressMessages(llply(c(20, 40, 60, 80, 100, 140), function(x) simbaltrees(x, metric="colless", numtrees=50, cutlow=-0.5, cuthigh=0.5), .parallel=TRUE))
out_larger_larger <- suppressMessages(llply(c(20, 40, 60, 80, 100, 140, 200, 300), function(x) simbaltrees(x, metric="colless", numtrees=50, cutlow=-0.5, cuthigh=0.5), .parallel=TRUE))

# Get traits
traitsoutdf <- ldply(out[[1]])

# Get network results
outdf <- ldply(out[[2]])
str(outdf)
head(outdf)
outdf_melt <- melt(outdf, id.vars=c(1,5,6))

outdf_melt_ <- ddply(outdf_melt, .(type, model, numsp, variable), summarise, 
                     mean = mean(value, na.rm=T),
                     se = sd(value, na.rm=T)/sqrt(na.omit(length(value)))
)
head(outdf_melt_)
splitted <- strsplit(outdf_melt_$model, " ")
repeater <- function(x) if(length(x)==1){rep(x,2)} else{x}
toadd <- ldply(splitted, repeater) # split names
outdf_melt_2 <- data.frame(outdf_melt_, toadd)
head(outdf_melt_2)

##### Visualize stuff
ggplot(outdf_melt_2[outdf_melt_2$variable=="connectance",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.5) +
  theme_bw() +
  facet_grid(V1 ~ V2, scales="free")

ggplot(outdf_melt_2[outdf_melt_2$variable=="nestedness",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.5) +
  theme_bw() +
  facet_grid(V1 ~ V2, scales="free")

# ggplot(outdf_melt_2[outdf_melt_2$variable=="links.per.species",], aes(numsp, mean, colour = type)) +
#   geom_point(size=4, alpha=0.5) +
#   theme_bw() +
#   facet_grid(V1 ~ V2, scales="free")


### 
outdf_larger <- ldply(out_larger)
str(outdf_larger)
head(outdf_larger)
outdf_larger_melt <- melt(outdf_larger, id.vars=c(1,7:9))
head(outdf_larger_melt)

outdf_larger_melt_ <- ddply(outdf_larger_melt, .(type, traitevol, netmethod, numsp, variable), summarise, 
                     mean = mean(value, na.rm=T),
                     se = sd(value, na.rm=T)/sqrt(na.omit(length(value)))
)

dfplot <- outdf_larger_melt_
ggplot(dfplot[dfplot$variable=="connectance",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.7) +
  theme_bw() +
  facet_grid(traitevol ~ netmethod, scales="free")
ggplot(dfplot[dfplot$variable=="nestedness",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.5) +
  theme_bw() +
  facet_grid(traitevol ~ netmethod, scales="free")
ggplot(dfplot[dfplot$variable=="web.asymmetry",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.7) +
  theme_bw() +
  facet_grid(traitevol ~ netmethod, scales="free")
ggplot(dfplot[dfplot$variable=="togetherness",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.7) +
  theme_bw() +
  facet_grid(traitevol ~ netmethod, scales="free")




### 
outdf_larger_larger <- ldply(out_larger_larger)
str(outdf_larger_larger)
head(outdf_larger_larger)
outdf_larger_larger_melt <- melt(outdf_larger_larger, id.vars=c(1,7:9))
head(outdf_larger_larger_melt)

outdf_larger_larger_melt_ <- ddply(outdf_larger_larger_melt, .(type, traitevol, netmethod, numsp, variable), summarise, 
    mean = mean(value, na.rm=T),
    se = sd(value, na.rm=T)/sqrt(na.omit(length(value))),
    n_trees = length(na.omit(value))
)
outdf_larger_larger_melt_$conf_int <- outdf_larger_larger_melt_$se*1.96
head(outdf_larger_larger_melt_)

dfplot <- outdf_larger_larger_melt_

limits <- aes(ymax = mean + conf_int, ymin = mean - conf_int)
dodge <- position_dodge(w=0.3)

ggplot(dfplot[dfplot$variable=="connectance",], aes(numsp, mean, shape = type, colour = type)) +
  geom_point(size=4, position=dodge) +
  geom_errorbar(limits, width=0.2, position=dodge) +
  theme_bw() +
  facet_grid(traitevol ~ netmethod, scales="free")
ggplot(dfplot[dfplot$variable=="nestedness",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.5) +
  theme_bw() +
  facet_grid(traitevol ~ netmethod, scales="free")
ggplot(dfplot[dfplot$variable=="web.asymmetry",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.7) +
  theme_bw() +
  facet_grid(traitevol ~ netmethod, scales="free")
ggplot(dfplot[dfplot$variable=="togetherness",], aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.7) +
  theme_bw() +
  facet_grid(traitevol ~ netmethod, scales="free")