library(ape)
library(phytools)
library(phangorn)
#install.packages("ips")
library(ips)

args = commandArgs(trailingOnly=TRUE)

trees.filename = args[1]
output.filename = args[2]

trees <- read.tree(trees.filename)

#collapse nodes
collapsed_tree_list <- list()
for (i in 1:length(trees)){
  tree = trees[[i]]
  tree_collapsed <- collapseUnsupportedEdges(tree, value = "node.label", 50 )
  collapsed_tree_list[[i]] <- tree_collapsed
}

class(collapsed_tree_list) <- "multiPhylo"

write.tree(collapsed_tree_list, file = paste0(output.filename,".tree"))
