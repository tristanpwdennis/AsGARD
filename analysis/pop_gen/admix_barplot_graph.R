library(tidyverse)
library(data.table)
library(starmie)
library(ape)
library(data.table)
library(phytools)
library(sf)
library(RColorBrewer)

# Load metadata and all admix data, find best Ks
df_samples <- read.csv('~/Projects/AsGARD/metadata/cease_combinedmetadata_noqc.20250212.csv')
df_samples$id <- seq(0, nrow(df_samples)-1)
#load admix path
admixpath <- '/Users/dennistpw/Projects/AsGARD/data/admixture_20250106/123/'

# log - likelihoods
LogLikelihood <- c(-21037075,-19628214,-18566335,-18201719,-18017098,-18052500,-17698206,-17425416.3,-17324632.0,-17112923, -16963538.0)
K <- seq(2, length(LogLikelihood)+1)
plot(x=K,y=LogLikelihood,type = "b", xaxt = "n")  # Log lik plateau around 6?, suggest optimal K is around there (eg 4 if sparse, 5 if comprehensive) as increases afterwards are quite marginal, we will plot only the first 10 then
axis(1, at = 1:12)
# Get haplonet output
filelist <- list.files('~/Projects/AsGARD/data/haplonet_20250211/admixture', full.names = TRUE)

#have to loop to get proper ordering
qfiles = list()
fileseq = seq(2, 7)
for (f in seq_along(fileseq)){
  qfiles[[f]] <- read.table(paste0('~/Projects/AsGARD/data/haplonet_20250211/admixture/CM023248.admixture.k',fileseq[[f]],'.q'))
}

# function for plotting admixture proportions 
plot_admix <- function(q, K){
  #get max index 
  find_max_index <- function(row) {
    return(which.max(row[1:K]))
  }

  # prepare admix df for k=5
admix_df <- q
admix_df$id <- seq(0, nrow(admix_df)-1)
admix_df$max_prob <- apply(admix_df, 1, find_max_index)
admix_df <- cbind(admix_df, df_samples)
admix_df <- admix_df %>% pivot_longer(cols=1:K, names_to = 'cluster', values_to = 'admixture_proportion')
admix_df$pop_code <- factor(admix_df$pop_code, levels=c("SAE", "SAR", "IRH","IRS", "APA", "INB","INM", "DJI","ETW","ETB","ETS","YEM","SUD"))
admixplot <- admix_df %>% mutate(id = forcats::fct_reorder(as.factor(id), as.numeric(max_prob))) %>% 
  ggplot(., aes(factor(id), admixture_proportion, fill = factor(cluster), color = NA)) +
  geom_col(color = NA, width = 1)+  # Set width to 1
  facet_grid(~pop_code, switch = "x", scales = "free", space = "free") +
  theme_minimal() +
  labs(y=K)+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(mult = 0.01))+
  scale_fill_brewer(palette = "Spectral")+
  theme(
    axis.text.x = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x.bottom = element_text(angle=90)
  ) 
return(admixplot)
}

 # hacky shit here as the order of files isn't in numeric - annoyingly

# Loop over sequence and plot all K's
admixplotlist <- list()
for (f in seq_along(fileseq)){
  admixplot = plot_admix(qfiles[[f]], fileseq[[f]])
  admixplotlist[[f]] <- admixplot
}

# Remove strip text bottom from all plots except the bottom
for (f in seq(1:6)) {
  admixplotlist[[f]] <- admixplotlist[[f]]+theme(strip.text.x.bottom = element_blank())
}

#remove panel labs from all except bottom plot

admix_multiK <- cowplot::plot_grid(plotlist = admixplotlist, ncol = 1, align = 'v')

ggsave(filename = '~/Projects/AsGARD/figures/admix_multiK.svg',width=10, height=6,)

#----------------------------------------------------------------------#
#Plotting Admix graph              #
#----------------------------------------------------------------------#

#define tip labels and colours
pop_cols <- c('#ff7f00','#507d2a','#007272','#33a02c','#a6cee3','#96172e','#f03e5e','#c57fc9','#c27a88','#6a3d9a','#cab2d6','#fccf86','#CC7722')

tipcols <- c('#007272','#33a02c','#fccf86','#CC7722','#a6cee3','#507d2a','#f03e5e','#96172e','#ff7f00','#c27a88','#c57fc9','#6a3d9a','#cab2d6')

#start by loading all trees
bootstrap_trees <- read.tree("~/Projects/AsGARD/data/treeMix_20250212/alltrees.treeout")
consensus_tree_bl <- consensus.edges(bootstrap_trees, p = 0.9) # Generate consensus tree
plot(consensus_tree_bl, type="unrooted")
node.labels <- round(prop.clades(consensus_tree_bl, bootstrap_trees) / 10) # Get bootstrap proportions
#write.tree(consensus_tree_bl, file = "~/Projects/AsGARD/data/treemix_20250106/consensus_unrooted.tre") # wo and root in figtree as this is a massive pain in r (node connecting APA with IRN/SAE)
#consensus_rooted <- read.nexus("~/Projects/AsGARD/data/treemix_20250106/consensus_rooted.tre")
svg(filename=paste0("~/Projects/AsGARD/figures/TreeMixOutput.svg",Sys.time(),".svg"), width = 6,height = 8) #Plot
ape::plot.phylo(consensus_tree_bl, tip.color =rev(tipcols), cex = 1, type = "unrooted")  # Plot with tip colors
nodelabels(node.labels, adj = c(1, -0.5), frame = "n", cex = 0.8, col = "black")
add.scale.bar(length=0.001)
dev.off()

