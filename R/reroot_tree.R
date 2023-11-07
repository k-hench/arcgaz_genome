library(phytools)
options(scipen = 6)
args <- commandArgs(trailingOnly = TRUE)

tree <- ape::read.tree(file = args[[1]])
tree_rerooted <- reroot(tree, node.number = as.numeric(args[[2]]))

tree_rotated <- tree_rerooted |>

ape::write.tree(phy = purrr::reduce(c(12, 18, 19, 20, 13, 14, 21, 17),
                                    ape::rotate,
                                    .init = tree_rerooted),
                file = args[[3]], digits = 6)