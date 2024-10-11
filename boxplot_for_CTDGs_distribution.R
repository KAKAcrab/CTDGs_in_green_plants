library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)
library(viridisLite)
library(viridis)
library(hrbrthemes)
library(ggbeeswarm)
library(ggpubr)
##boxplot for distribution of CTDGs in chromosome of each species

box_data <- read_excel("~/R  project/207spe_ctdg/sorted_distribution_of_CTDGs_in_chromosome.xlsx")
box_data$species <- factor(box_data$species,levels = rev(unique(box_data$species)))
box_data <- box_data %>% distinct()
head(box_data)
light_colors <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462",
                  "#b3de69","#a6cee3","#1f78b4","#fb9a99","#6a3d9a","#b15928","#e31a1c")
p<-ggplot(box_data, aes(x=count,y=species)) + 
  geom_boxplot(cex = .1,
               outlier.size =.5,
               # corral = "gutter",
               # corral.width = 0.5,
               # method = "center",
               # size = .5,
               color="black",
               show.legend=T) +
  # geom_violin(trim = FALSE,
  #             alpha = 1,
  # 
  #             width = 1,
  #             #adjust=0.5,
  #             color="black",
  #             size = 0.2,
  #             show.legend = FALSE)+
  # geom_jitter(shape=21,position = position_jitter(width = 0.2),color="black")+
  # ylim(0,700)+
  # scale_fill_manual(values=light_colors) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x=element_line(color = "black"),
        axis.text.x=element_text(angle = 0,hjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 2.5,color = "black"),
        strip.text = element_text(size = 8),#调整每个分组标题的字体大小
        legend.position = "top"
  )+
  # legend.position = "none") +
  guides(fill=guide_legend(),color="none")+
  # coord_flip() +
  labs(
    x = "",
    # x="group",
    y = ""
    # y= "Number of genes"
    #title = "Distribution of CTDGs across the 7 pear species"
  );p
  # facet_wrap(~group, nrow = 3)+
 
setwd("~/Desktop/project/CTDG_of_207_spe/draft_fig/")
# ggsave("boxplot_of_taxon_5*2.pdf", plot = p, width = 5, height = 2, units = "in")