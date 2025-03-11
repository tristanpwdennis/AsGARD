library(tidyverse)
library(data.table)
library(RColorBrewer)
library(viridis)


# Load data
cnv_iran <- fread('/Users/dennistpw/Projects/AsGARD/data/cnv/modal_copy_number.csv', header=T)
sample_code_iran <- fread('/Users/dennistpw/Projects/AsGARD/data/cnv/iransamples.csv', header=TRUE)

cnv_noiran <-  fread('~/Projects/cease/cease_wp1c/analysis_data/cnv_calling/modal_copy_number_all.tsv', header=T)
cnv_iran <- cbind(pop_code=sample_code_iran$pop_code, cnv_iran)
cnv_iran <- cbind(analysis_pop = 'iran', cnv_iran)
cnv <- rbind(cnv_noiran, cnv_iran)
cnv_genes <- fread('~/Projects/AsGARD/data/cnv_gene_locs.csv', header=T)
#Select genes of interest from modal CNV table, pivot to long
geneids <- cnv_genes[,c('AnsteID')]
cnv_wide <- subset(cnv, select = names(cnv) %in% geneids$AnsteID)
cnv_wide <- cbind(cnv$sample_id,cnv$pop_code, cnv_wide)
cnv_long <- pivot_longer(cnv_wide, cols = 3:31)
colnames(cnv_long) <- c('sample_id','pop_code','AnsteID','copynumber')
cnv_long <- left_join(cnv_long, cnv_genes, by='AnsteID')

# Order populations for plotting
cnv_long$analysis_pop = factor(cnv_long$pop_code , levels=c("SAE","SAR","IRH",'IRS','APA',"INB","INM",'DJI','ETW','ETB','ETS','SUD','YEM'), labels=c("SAE","SAR","IRH",'IRS','APA',"INB","INM",'DJI','ETW','ETB','ETS','SUD','YEM')) 
cnv_long<-cnv_long[ order(cnv_long$AnsteChrom, cnv_long$AnsteStart), ]
cnv_long$gene_id <- as.factor(cnv_long$AnsteID)

# Plot

cnv_long_ordered <- cnv_long %>%
  group_by(AnsteID, AnsteFamily) %>%
  arrange(AnsteStart, .by_group = TRUE) %>%
  ungroup()
cnv_heatmap <- cnv_long_ordered %>% #filter(AnsteFamily != 'cyp9k1' & AnsteFamily != "coejh" & AnsteFamily != 'ace1') %>% 
  ggplot(aes(x=as.factor(sample_id),y = reorder(GeneName, AnsteStart), fill=copynumber))+
  geom_tile()+
  facet_grid(rows=vars(AnsteFamily), cols =vars(analysis_pop),scales="free",space="free")+
  theme_minimal()+
  scale_fill_distiller(palette = 'Blues', direction = 1)+
  labs(fill="Est. Copy Number", y='Gene',x='Individual')+
  theme(
    axis.ticks.x = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    #axis.text.y = element_blank(),
    strip.text.x = element_text(size=10,face="bold", angle=90), #
    #axis.text.y=element_text(size=15,face="bold"),
    strip.text.y = element_text(size=12,face="bold",angle = 90),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.spacing = element_blank(),
    plot.margin = margin(r = 0.05, l = 0.05))

cnv_heatmap
ggsave('~/Projects/AsGARD/figures/cnv_heatmap.svg', plot = cnv_heatmap, width = 9, height = 5)
ggsave('~/Projects/AsGARD/figures/cnv_heatmap.png', plot = cnv_heatmap, width = 9, height = 5, bg = NULL)

## Now generate SNP heatmap
snp_data <- fread('~/Projects/AsGARD/data/irgenes_snps.csv')
snp_data <- snp_data %>% filter(gene %in% c('Ace1','Vgsc','Rdl'))
snp_data_long <- pivot_longer(snp_data, cols = 2:13) 
snp_data_long$name <- factor(snp_data_long$name,levels=c("SAE","SAR","IRH",'IRS','APA',"INB","INM",'DJI','ETW','ETB','ETS','SUD','YEM'), labels=c("SAE","SAR","IRH",'IRS','APA',"INB","INM",'DJI','ETW','ETB','ETS','SUD','YEM')) 
snp_data_long$value <- as.numeric(snp_data_long$value)


snp_heatmap <- snp_data_long %>% ggplot(aes(x=name,y=effect, fill=value))+
  facet_grid(rows=vars(gene), scales = 'free_y', space='free_y')+
  geom_tile()+
  geom_text(aes(label = value), colour='darkgrey')+
  scale_fill_gradient(high = "#132B43", low = "white")+
  #scale_color_gradient(low = "#132B43", high = "white")+
  theme_classic()+
  theme(
    panel.spacing = unit(0.5, "lines"), 
    strip.text.y = element_text(angle = 270), 
    legend.position = 'none'
  ) +
  labs(x = "Cohort", y = "Variant Effect", fill = "Frequency")
ggsave(filename='~/Projects/AsGARD/figures/snp_heatmap.png',width=6, height=9)
