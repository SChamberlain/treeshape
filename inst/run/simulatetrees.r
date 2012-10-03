################## Simulate trees using function simbal
install.packages(c("plyr","ggplot2","ape","apTreeshape","bipartite","reshape2"))
library(plyr); library(ape); library(apTreeshape); library(bipartite); library(ggplot2); library(reshape2)
source("fastBM.r"); source("nodeHeights.r"); source("pbtree.r")

### plant trees
foo <- function(tips) {
  
  trees_colless_plants <- simbal(t=tips, "colless", 1000, cutlow = -0.8, cuthigh = 0.8)
  trees_colless_plants_bal <- compact(lapply(trees_colless_plants, function(x) x$colless_bal)) # get the balanced trees
  trees_colless_plants_unbal <- compact(lapply(trees_colless_plants, function(x) x$colless_unbal)) # get the unbalanced trees
  
  ### animal trees
  trees_colless_anim <- simbal(t=tips, "colless", 1000, cutlow = -0.8, cuthigh = 0.8)
  trees_colless_anim_bal <- compact(lapply(trees_colless_anim, function(x) x$colless_bal)) # get the balanced trees
  trees_colless_anim_unbal <- compact(lapply(trees_colless_anim, function(x) x$colless_unbal)) # get the unbalanced trees
  
  ################## Simulate traits on each tree
  # Plants
  t1_col_plants_bal <- llply(trees_colless_plants_bal, fastBM, a = 10, bounds=c(0,100))
  t1_col_plants_unbal <- llply(trees_colless_plants_unbal, fastBM, a = 10, bounds=c(0,100))
  t2_col_plants_bal <- llply(trees_colless_plants_bal, fastBM, a = 10, bounds=c(0,100))
  t2_col_plants_unbal <- llply(trees_colless_plants_unbal, fastBM, a = 10, bounds=c(0,100))
  abd_col_plants_bal <- round(rlnorm(length(trees_colless_plants_bal), meanlog=5, 1), 0)
  abd_col_plants_unbal <- round(rlnorm(length(trees_colless_plants_unbal), meanlog=5, 1), 0)
  
  # Animals
  t1_col_anim_bal <- llply(trees_colless_anim_bal, fastBM, a = 10, bounds=c(0,100))
  t1_col_anim_unbal <- llply(trees_colless_anim_unbal, fastBM, a = 10, bounds=c(0,100))
  t2_col_anim_bal <- llply(trees_colless_anim_bal, fastBM, a = 10, bounds=c(0,100))
  t2_col_anim_unbal <- llply(trees_colless_anim_unbal, fastBM, a = 10, bounds=c(0,100))
  abd_col_anim_bal <- round(rlnorm(length(trees_colless_anim_bal), meanlog=5, 1), 0)
  abd_col_anim_unbal <- round(rlnorm(length(trees_colless_anim_unbal), meanlog=5, 1), 0)
  
  
  ################## Get plant-animal phylogenetic tree pairs
  # Balanced trees
  all_bal <- list(trees_colless_plants_bal, trees_colless_anim_bal)
  smaller_bal <- which.min(sapply(all_bal, length))
  larger_bal <- which.max(sapply(all_bal, length))
  tree_pairs_bal <- list(all_bal[[smaller_bal]], all_bal[[larger_bal]][1:length(all_bal[[smaller_bal]])])
  
  # Unbalanced trees
  all_unbal <- list(trees_colless_plants_unbal, trees_colless_anim_unbal)
  smaller_unbal <- which.min(sapply(all_unbal, length))
  larger_unbal <- which.max(sapply(all_unbal, length))
  tree_pairs_unbal <- list(all_unbal[[smaller_unbal]], all_unbal[[larger_unbal]][1:length(all_unbal[[smaller_unbal]])])
  
  # Pairs of traits from balanced trees
  all_t1_bal <- list(t1_col_plants_bal[1:length(tree_pairs_bal[[1]])], 
                     t1_col_anim_bal[1:length(tree_pairs_bal[[1]])])
  all_t2_bal <- list(t2_col_plants_bal[1:length(tree_pairs_bal[[1]])], 
                     t2_col_anim_bal[1:length(tree_pairs_bal[[1]])])
  all_abd_bal <- list(abd_col_plants_bal[1:length(tree_pairs_bal[[1]])], 
                      abd_col_anim_bal[1:length(tree_pairs_bal[[1]])])
  
  # Pairs of traits form unbalanced trees
  all_t1_unbal <- list(t1_col_plants_unbal[1:length(tree_pairs_unbal[[1]])], 
                       t1_col_anim_unbal[1:length(tree_pairs_unbal[[1]])])
  all_t2_unbal <- list(t2_col_plants_unbal[1:length(tree_pairs_unbal[[1]])], 
                       t2_col_anim_unbal[1:length(tree_pairs_unbal[[1]])])
  all_abd_unbal <- list(abd_col_plants_unbal[1:length(tree_pairs_unbal[[1]])], 
                        abd_col_anim_unbal[1:length(tree_pairs_unbal[[1]])])
  
  ################## Simulate networks on each tree
  # Simulate completely random networks, regardless of traits, etc.
  mats_rand_bal <- sim_rand_nets(tree_pairs_bal)
  mats_rand_unbal <- sim_rand_nets(tree_pairs_unbal)
  
  # Simulate networks, with interactions propoprtional to trait matching
  # ratio
  mats_traits1_bal_ratio <- sim_traits_nets(all_t1_bal, method = "r", value = 1.5)
  mats_traits1_unbal_ratio <- sim_traits_nets(all_t1_unbal, method = "r", value = 1.5)
  mats_traits2_bal_ratio <- sim_traits_nets(all_t2_bal, method = "r", value = 1.5)
  mats_traits2_unbal_ratio <- sim_traits_nets(all_t2_unbal, method = "r", value = 1.5)
  mats_traitsabd_bal_ratio <- sim_traits_nets(all_abd_bal, method = "r", value = 1.5)
  mats_traitsabd_unbal_ratio <- sim_traits_nets(all_abd_unbal, method = "r", value = 1.5)
  
  # complementarity
  mats_traits1_bal_comp <- sim_traits_nets(all_t1_bal, method = "c", value = 2)
  mats_traits1_unbal_comp <- sim_traits_nets(all_t1_unbal, method = "c", value = 2)
  mats_traits2_bal_comp <- sim_traits_nets(all_t2_bal, method = "c", value = 2)
  mats_traits2_unbal_comp <- sim_traits_nets(all_t2_unbal, method = "c", value = 2)
  mats_traitsabd_bal_comp <- sim_traits_nets(all_abd_bal, method = "c", value = 2)
  mats_traitsabd_unbal_comp <- sim_traits_nets(all_abd_unbal, method = "c", value = 2)
  
  # barrier (no need to give value parameter)
  mats_traits1_bal_barr <- sim_traits_nets(all_t1_bal, method = 'b')
  mats_traits1_unbal_barr <- sim_traits_nets(all_t1_unbal, method = 'b')
  mats_traits2_bal_barr <- sim_traits_nets(all_t2_bal, method = 'b')
  mats_traits2_unbal_barr <- sim_traits_nets(all_t2_unbal, method = 'b')
  mats_traitsabd_bal_barr <- sim_traits_nets(all_abd_bal, method = 'b')
  mats_traitsabd_unbal_barr <- sim_traits_nets(all_abd_unbal, method = 'b')
  
  
  ################## Calculate network metrics on matrices
  df_rand <- getnetmets(mats_rand_bal, mats_rand_unbal) # random networks
  
  df_traits1_ratio <- getnetmets(mats_traits1_bal_ratio, mats_traits1_unbal_ratio) # ratio traits
  df_traits1_comp <- getnetmets(mats_traits1_bal_comp, mats_traits1_unbal_comp) # complementarity traits
  df_traits1_barr <- getnetmets(mats_traits1_bal_barr, mats_traits1_unbal_barr) # barrier traits
  
  df_traits2_ratio <- getnetmets(mats_traits2_bal_ratio, mats_traits2_unbal_ratio) # ratio traits
  df_traits2_comp <- getnetmets(mats_traits2_bal_comp, mats_traits2_unbal_comp) # complementarity traits
  df_traits2_barr <- getnetmets(mats_traits2_bal_barr, mats_traits2_unbal_barr) # barrier traits
  
  alldat <- rbind(df_rand, df_traits1_ratio, df_traits2_ratio,
                  df_traits1_comp, df_traits2_comp,
                  df_traits1_barr, df_traits2_barr)
  alldat$model <- c( rep("Random",nrow(df_rand)), 
                     rep("Ratio 1",nrow(df_traits1_ratio)),
                     rep("Ratio 2",nrow(df_traits2_ratio)),
                     rep("Complementarity 1",nrow(df_traits1_comp)),
                     rep("Complementarity 2",nrow(df_traits2_comp)),
                     rep("Barrier 1",nrow(df_traits1_barr)),
                     rep("Barrier 2",nrow(df_traits2_barr))
  )
  alldat$numsp <- rep(tips, nrow(alldat))
  alldat
}

out <- suppressMessages(llply(c(20, 40, 60, 80, 100, 140, 160), foo, .progress="text"))
outdf <- ldply(out)
outdf_melt <- melt(outdf, id.vars=c(1,3,4))

outdf_melt_ <- ddply(outdf_melt, .(type, model, numsp, variable), summarise, 
                     mean = mean(value, na.rm=T),
                     se = sd(value, na.rm=T)/sqrt(na.omit(length(value)))
)

##### Visualize stuff
ggplot(outdf_melt_, aes(numsp, mean, colour = type)) +
  geom_point(size=4, alpha=0.5) +
  theme_bw() +
  facet_wrap(~ model)