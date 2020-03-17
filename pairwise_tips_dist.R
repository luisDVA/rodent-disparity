library(dispRity) # Measuring Disparity
library(Rphylopars) # Phylogenetic Comparative Tools for Missing Data and Within-Species Variation

set.seed(12)
rand_tree <- rtree(6) 
rand_tree <- force.ultrametric(rand_tree, method=c("extend"))
#rand_tree$node.label <- paste0("n", 5:7)
plot(rand_tree)
nodelabels()
dummy_matrix <- matrix(rnorm(24,2), 6, 4)
dummy_matrix
row.names(dummy_matrix) <- rand_tree$tip.label
dummy_matrix
node_vals <- anc.recon(dummy_matrix, rand_tree, vars = FALSE, CI = FALSE)
tipsnodesmat <- rbind(dummy_matrix,node_vals)
tipsnodesmat
#rownames(tipsnodesmat) <- c(rand_tree$tip.label, rand_tree$node.label)
anc_dists <- ancestral.dist(matrix = tipsnodesmat, tree = rand_tree,full = FALSE)
anc_dists



pairwise_tip_distance <- function(tree, tips) {
  anc_dists <- ancestral.dist(matrix = tipsnodesmat, tree = rand_tree,full = FALSE)
  tips_mrca <- getMRCA(tree, tips)
  nodenums <- seq((length(tree$tip.label) + 1), ((length(tree$tip.label)) + tree$Nnode), 1)
  nodedescs <- phylobase::descendants(as(tree, "phylo4"), tips_mrca, "all")
  int_node <- unname(nodedescs[which(is.na(names(nodedescs)))])
  distvecsum <- c(tips, int_node)
    sum(anc_dists[distvecsum])
}

library(purrr) # Functional Programming Tools
library(RcppAlgos) # High Performance Tools for Combinatorics and Computational Mathematics
library(phytools) # Phylogenetic Tools for Comparative Biology (and Other Things)

tipcombos <-  as.data.frame(comboGeneral(rand_tree$tip.label, m=2),stringsAsFactors = FALSE)
tipcomboslist <- purrr::map(transpose(tipcombos),paste)

pw_distances <- map_dbl(tipcomboslist,~pairwise_tip_distance(tree = rand_tree,tips = .x))
pw_distances

# comparing distances
anc_dists_tips <- tibble::enframe(anc_dists)
anc_dists_tips <- anc_dists_tips[anc_dists_tips$name %in% rand_tree$tip.label,]
anc_dists_tips  
sister_dists <- cbind(tipcombos,pw_distances)
sister_dists



library(ggplot2)  
library(dplyr)
ancDF <- anc_dists_tips %>% rename(tip=1,dist=2) %>% mutate(dist_type="ancestral")
sistDF <- sister_dists %>% tidyr::unite("tip",V1:V2) %>% rename(dist=pw_distances) %>% mutate(dist_type="tip")
dispDistrib <- bind_rows(ancDF,sistDF)
  
ggplot(dispDistrib)+
  geom_boxplot(aes(x=dist_type,y=dist))+
  geom_jitter(aes(x=dist_type,y=dist))

# data.frame(tip=rand_tree$tip.label,
#            sister= simplify(flatten(purrr::map(rand_tree$tip.label,~getSisters(rand_tree,.x, mode = "label")))),
#            stringsAsFactors = FALSE)






