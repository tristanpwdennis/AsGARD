library(tidyverse)
library(ape)
library(ggtree)
library(treeio)

# col mapping dict
tipcols <- c('ETB'='#007272',
             'ETS'= '#33a02c',
             'SUD' = '#fccf86',
             'YEM'= '#CC7722',
             'ETW' = '#a6cee3',
             "DJI"= '#507d2a',
             "INM" = '#f03e5e',
             "INB" = '#96172e',
             "APA" = '#ff7f00',
             "IRH" = '#c57fc9',
             "IRS" = '#c27a88',
             "SAE" = '#6a3d9a',
             "SAR" = '#cab2d6',
             "MACULATUS" = 'lightgrey')


# Read while mito tree
tree <- read.tree('~/Projects/AsGARD/data/mitochondrial_phylogeny/consensus_sequences.fasta.treefile') 
tree <- root(tree, outgroup = "KT382822.1", resolve.root = TRUE) #root on macultaus
tree$edge.length[tree$edge[,1] == Ntip(tree) + 1] <- tree$edge.length[tree$edge[,1] == Ntip(tree) + 1] / 100  # Shorten root

tree_ladderized <- ladderize(tree)  # tidy with ladderise

#load annotations
annotations <- read.csv('~/Projects/AsGARD/data/mitochondrial_phylogeny/tree_annotations.csv')
all(tree_ladderized$tip.label %in% annotations$label)  # do content match?
# Reorder annotations based on tree tip labels
annotations <- annotations[match(tree_ladderized$tip.label, annotations$label), ]
head(annotations)

#Plot
svg(filename=paste0("~/Projects/AsGARD/figures/mito_whole_tree.svg",Sys.time(),".svg"), width = 6,height = 8) #Plot
plot(tree_ladderized, show.tip.label = FALSE)
tiplabels(pch = 19, col = tipcols[annotations$cohort], cex = 0.5)
add.scale.bar()
dev.off()


coitree <- read.tree('~/Projects/AsGARD/data/mitochondrial_phylogeny/coi.fasta.treefile') 
#coitree <- root(coitree, outgroup = "MACULATUS_PUB", resolve.root = TRUE) #root on macultaus
coitree$edge.length[coitree$edge[,1] == Ntip(coitree) + 1] <- coitree$edge.length[coitree$edge[,1] == Ntip(coitree) + 1] / 100  # Shorten root
coitree <- ladderize(coitree)
plot(coitree)

#load coi anns
coi_ann <- read.csv('~/Projects/AsGARD/data/mitochondrial_phylogeny/coi_tree_annotations.csv')
all(coitree$tip.label %in% coi_ann$label)  # do content match?
coi_ann <- coi_ann[match(coitree$tip.label, coi_ann$label), ]

# Sort tip colours
coitipcols <- tipcols[coi_ann$cohort]
coitipcols[is.na(coitipcols)] <- "grey"

modified_labels <- ifelse(grepl('PUB', coitree$tip.label), coitree$tip.label, "")
coitree$tip.label <- modified_labels

svg(filename=paste0("~/Projects/AsGARD/figures/coi_whole_tree.svg",Sys.time(),".svg"), width = 6,height = 8) #Plot
plot(coitree, show.tip.label = TRUE,cex.tip = 0.6 )
tiplabels(pch = 19, col = coitipcols, cex = 0.5)
add.scale.bar()
dev.off()

