library(phytools)
options(scipen = 6)
args <- commandArgs(trailingOnly = TRUE)

tree <- ape::read.tree(file = args[[1]])
tree_rerooted <- reroot(tree, node.number = as.numeric(args[[2]]))

tree_rotated <- tree_rerooted |>
    ape::rotate(12) |>
    ape::rotate(18) |>
    ape::rotate(19) |>
    ape::rotate(20) |>
    ape::rotate(13) |>
    ape::rotate(14) |>
    ape::rotate(21) |>
    ape::rotate(17)

ape::write.tree(phy = tree_rotated,
                file = args[[3]], digits = 6)