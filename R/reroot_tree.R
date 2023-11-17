library(phytools)
options(scipen = 6)
args <- commandArgs(trailingOnly = TRUE)

tree <- ape::read.tree(file = args[[1]])
tree_rerooted <- reroot(tree, node.number = as.numeric(args[[2]]))

ape::write.tree(phy = purrr::reduce(c(12,18,21,17,19,20,21,14,15),
                                    ape::rotate,
                                    .init = tree_rerooted),
                file = args[[3]], digits = 6)
