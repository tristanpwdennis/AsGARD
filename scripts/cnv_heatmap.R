library(tidyverse)
library(data.table)
library(RColorBrewer)
library(viridis)

# Load data
cnv <- fread('~/Projects/cease/cease_wp1c/analysis_data/cnv_calling/modal_copy_number_all.tsv', header=T)
cnv_genes <- fread('~/Projects/AsGARD/data/cnv_gene_locs.csv', header=T)

#Select genes of interest from modal CNV table, pivot to long
geneids <- cnv_genes[,c('AnsteID')]
cnv_wide <- subset(cnv, select = names(cnv) %in% geneids$AnsteID)
cnv_wide <- cbind(cnv$sample_id,cnv$pop_code, cnv_wide)
cnv_long <- pivot_longer(cnv_wide, cols = 3:31)
colnames(cnv_long) <- c('sample_id','pop_code','AnsteID','copynumber')
cnv_long <- left_join(cnv_long, cnv_genes, by='AnsteID')

# Order populations for plotting
cnv_long$analysis_pop = factor(cnv_long$pop_code , levels=c("SAE","SAR","IRN","INB","INM","APA",'DJI','ETW','ETB','ETS','SUD','YEM'), labels=c("SAE","SAR","IRN","INB","INM","APA",'DJI','ETW','ETB','ETS','SUD','YEM')) 
cnv_long<-cnv_long[ order(cnv_long$AnsteChrom, cnv_long$AnsteStart), ]
cnv_long$gene_id <- as.factor(cnv_long$AnsteID)

# Plot

cnv_long_ordered <- cnv_long %>%
  group_by(AnsteID, AnsteFamily) %>%
  arrange(AnsteStart, .by_group = TRUE) %>%
  ungroup()
cnv_heatmap <- cnv_long_ordered %>% filter(AnsteFamily != 'cyp9k1' & AnsteFamily != "coejh" & AnsteFamily != 'ace1') %>% 
  ggplot(aes(x=as.factor(sample_id),y = reorder(GeneName, AnsteStart), fill=copynumber))+
  geom_tile()+
  facet_grid(rows=vars(AnsteFamily), cols =vars(analysis_pop),scales="free",space="free", switch = 'x')+
  theme_minimal()+
  scale_fill_distiller(palette = 'Blues', direction = 1)+
  labs(fill="Est. Copy Number", y='Gene',x='Individual')+
  theme(
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        #axis.text.y = element_blank(),
        strip.text.x.bottom = element_text(size=18,face="bold", angle=90), #THE LINE THAT IS EDITED
        #axis.text.y=element_text(size=15,face="bold"),
        strip.text.y = element_text(size=18,face="bold",angle = 0),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.spacing = element_blank(),
        plot.margin = margin(r = 0.05, l = 0.05))

cnv_heatmap
ggsave('~/Projects/AsGARD/figures/cnv_heatmap.svg', plot = cnv_heatmap, width = 20, height = 7, bg = NULL)
ggsave('~/Projects/AsGARD/figures/cnv_heatmap.png', plot = cnv_heatmap, width = 20, height = 7, bg = NULL)


