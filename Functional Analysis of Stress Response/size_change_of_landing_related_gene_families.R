library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)
library(viridisLite)
library(viridis)
library(hrbrthemes)
library(ggbeeswarm)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
################################################################################################################
################################################################################################################
desired_order_1 <- c(
  "Angiosperms","Gymnosperms","Ferns","Lycophytes","Bryophytes","Zygnematophyceae","Streptophyte_algae","Chlorophyta"
)
desired_order_2 <- c(
  "Lycophytes","Bryophytes","Zygnematophyceae","Streptophyte_algae","Chlorophyta"
)
my_comparisons_1 <- list(
  c("Lycophytes","Bryophytes"),c("Bryophytes","Zygnematophyceae"),
  c("Zygnematophyceae","Streptophyte_algae"),c("Streptophyte_algae","Chlorophyta"),
  c("Ferns","Lycophytes"),c("Gymnosperms","Ferns"),c("Angiosperms","Gymnosperms")
)
my_comparisons_2 <- list(
  c("Lycophytes","Bryophytes"),c("Bryophytes","Zygnematophyceae"),
  c("Zygnematophyceae","Streptophyte_algae"),c("Streptophyte_algae","Chlorophyta")
)
light_colors <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462",
                  "#b3de69","#a6cee3","#1f78b4","#fb9a99","#6a3d9a","#b15928","#e31a1c")

################################################################################################################
################################################################################################################
#各类群中CTDGs的数量
box_data_1 <- read_excel("~/Desktop/project/CTDG_of_207_spe/genes_involved_in_landing/number_of_clustered_genes_involved_in_landing.xlsx",
                       sheet="number_of_CTDGs_across_7_taxon")
box_data_1$Species <- factor(box_data_1$Species,levels =rev(desired_order_1))
box_data_1 <- box_data_1 %>% distinct()

p1 <-ggplot(box_data_1, aes(x=Species,y=value,color=Species)) + 
  geom_boxplot(cex = .2,
               outlier.size =.5,
               # corral = "gutter",
               # corral.width = 0.5,
               # method = "center",
               # size = .5,
               # color="black",
               show.legend=T) +
  # ylim(0,3000)+
  scale_fill_manual(values=light_colors) +
  stat_compare_means(
                     # aes(group=Species),
                     comparisons=my_comparisons_1,
                     method = "wilcox.test",
                     # method = "t.test",
                     exact=F,
                     label = "p.signif",
                     show.legend = F,
                     size=2
                     ) +
  theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.grid.major.x = element_line(color = "black",size = .2,linetype = "dashed"),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_text(angle = 0,hjust = 0.5,size = 7,color = "black"),
        axis.text.y = element_text(size = 7,color = "black"),
        axis.ticks.x = element_line(color = "black",size = .1),
        axis.ticks.y = element_line(color = "black",size = .1),
        panel.border = element_rect(color = "black", fill = NA, size = .3), # 增加边框
        strip.text = element_text(size = 8),#调整每个分组标题的字体大小
        legend.position = "top"
  )+
  # legend.position = "none") +
  guides(fill=guide_legend(),color="none")+
  coord_flip() +
  labs(
    x = "",
    # x="group",
    y = ""
    # y= "Number of genes"
    #title = "Distribution of CTDGs across the 7 pear species"
  ) +
  # facet_wrap(~group, nrow = 3)+
  scale_x_discrete(limits = rev(desired_order_1));p1
################################################################################################################
################################################################################################################
##登陆性状成簇基因的数量
box_data_2 <- read_excel("~/Desktop/project/CTDG_of_207_spe/genes_involved_in_landing/number_of_clustered_genes_involved_in_landing.xlsx",
                       sheet = "boxplot")
box_data_2 <- gather( box_data_2,group,value,-Species)
box_data_2$Species <- factor(box_data_2$Species,levels =rev(c("Lycophytes","Bryophytes","Zygnematophyceae",
                                                  "Streptophyte_algae","Chlorophyta")) )
box_data_2$group <- factor(box_data_2$group, levels = c("Abscisic acid","formation of plant cuticle",
                                                        "UV-B","cold tolerance","cell wall",
                                                        "heat tolerance","salt tolerance"))
# view(box_data_2)
p2 <-ggplot(box_data_2, aes(x=Species,y=value,color=Species)) + 
  geom_boxplot(cex = .2,
               outlier.size =.5,
               # corral = "gutter",
               # corral.width = 0.5,
               # method = "center",
               # size = .5,
               # color="black",
               show.legend=T) +
  # ylim(0,2000)+
  scale_fill_manual(values=light_colors) +
  stat_compare_means(
                     # aes(group=Species),
                     comparisons=my_comparisons_2,
                     method = "wilcox.test",
                     # method = "t.test",
                     exact=F,
                     label = "p.signif",
                     show.legend = F,
                     size=2
                     ) +
  theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.grid.major.x = element_line(color = "black",size = .2,linetype = "dashed"),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        axis.text.x=element_text(angle = 0,hjust = .5,vjust=.5,size = 7,color = "black"),
        axis.text.y = element_text(hjust = 1,vjust=.5,size = 7,color = "black"),
        axis.ticks.x = element_line(color = "black",size = .1),
        axis.ticks.y = element_line(color = "black",size = .1),
        strip.text = element_text(size = 8,color = "black",face = "bold",hjust = 0),#调整每个分组标题的字体大小
        strip.background = element_blank(), # 隐藏子图标题部分的灰色背景
        panel.border = element_rect(color = "black", fill = NA, size = .3), # 增加边框
        legend.position = "top"
  )+
  # legend.position = "none") +
  guides(fill=guide_legend(),color="none")+
  coord_flip() +
  labs(
    x = "",
    # x="group",
    y = ""
    # y= "Number of genes"
    #title = "Distribution of CTDGs across the 7 pear species"
  ) +
  facet_wrap(~group, nrow = 2)+
  scale_x_discrete(limits = rev(desired_order_2));p2
################################################################################################################
################################################################################################################
family_levels <- colnames(read_excel("~/Desktop/project/CTDG_of_207_spe/genes_involved_in_landing/number_of_clustered_genes_involved_in_landing.xlsx",
                                       sheet="genefamily",col_names = TRUE))[-1]
heatmap_data_1<- read_excel("~/Desktop/project/CTDG_of_207_spe/genes_involved_in_landing/number_of_clustered_genes_involved_in_landing.xlsx",
                            sheet="genefamily",col_names = TRUE)
head(heatmap_data_1)
heatmap_data_1 <- gather(heatmap_data_1,key = "family",value="Count",-"Species")
# heatmap_data_1 <- merge(heatmap_data_1,group_info,by="family",all=F)
heatmap_data_1$family <- factor(heatmap_data_1$family,levels =family_levels)
# heatmap_data_1$Count <- as.numeric(heatmap_data_1$Count)
# heatmap_data_1 <- as.data.frame(heatmap_data_1)
head(heatmap_data_1)
p3 <- ggplot(heatmap_data_1%>% filter(Count != 0), aes(y = family, x = Species)) +
  geom_point(aes(size = Count,fill=family), shape = 21, color = "black", alpha = 0.7,stroke=0.5) +
  scale_size_area(max_size = 11) +
  scale_fill_manual(values = light_colors)+
  theme_minimal() +
  # coord_flip() +  # 坐标轴翻转
  labs(x = "",
       y = "",
       size = "Count",
       tag = "") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0,hjust=0.5,size = 7,face = "bold",color = "black"),
        axis.text.y = element_text(angle = 0, size = 7,color = "black"),
        # axis.line.x = element_line(color = "black"), 
        axis.line.x = element_blank(),
        # axis.ticks.x.top = element_line(color = "black"), ##保留刻度线
        # axis.line.y=element_line(color = "black"),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        #legend.text = element_text(angle = 90),
        panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        # panel.border = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.2),
        axis.line = element_line(color = "black")
  )+
  guides(fill = F)+  #隐藏代表物种颜色的图例
  scale_y_discrete(limits = rev(family_levels),position = "left")+
  scale_x_discrete(position = "top");p3
################################################################################################################
################################################################################################################
setwd("~/Desktop/project/CTDG_of_207_spe/draft_fig/")
p <- plot_grid(p1,p3,p2, nrow = 3,rel_heights = c(8,10,12),rel_widths = c(10,10,10));p
ggsave("size_change_of_landing_related_gene_families_8.26*9.pdf", plot = p, width = 8.26, height = 9, units = "in")





















