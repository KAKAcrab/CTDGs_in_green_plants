library(pheatmap)
library(readxl)
library(tibble)
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(reshape2)
library(scales)
library(hrbrthemes)

colors_2 <-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
             "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
             "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")
##气泡图 
data_3 <- read_excel("~/R  project/207spe_ctdg/gene_family_counts.xlsx",sheet="mean")
data_3$LatinName <- factor(data_3$LatinName,levels = rev(unique(data_3$LatinName)))
data_3 <- gather(data_3,genefamily,count,-LatinName)
data_3 <- data_3%>%
  mutate(count = as.numeric(count))
# %>%
#   mutate(count = log2(count + 1))
# View(data_4)

p3 <- ggplot(data_3, aes(y = genefamily, x = LatinName)) +
  geom_blank() +
  # ggplot(data = subset(data_4, count > 0), aes(y = genefamily, x = LatinName)) +
  geom_point(data = subset(data_3, count > 0),
             aes(size = count,fill=LatinName),
             shape = 21, color = "black", alpha = 0.7,stroke=0.5) +
  scale_size_area(max_size = 8) +
  scale_fill_manual(values = rev(colors_2)) +  # 自定义颜色映射
  #scale_fill_manual(values = light_colors)+
  theme_minimal() +
  # coord_flip() +  # 坐标轴翻转
  labs(x = "",
       y = "",
       size = "Count",
       tag = "") +
  theme(axis.text.x = element_text(angle = 30, vjust = 1,hjust=1,size = 9,color = "black"),
        axis.text.y = element_text(angle = 0, size = 8,color = "black"),
        axis.line.x = element_blank(), ##去除坐标线
        axis.ticks.x.top = element_line(color = "black"), ##保留刻度线
        # axis.line.y=element_line(color = "black"),
        axis.line.y = element_blank(),
        legend.position = "right",
        # legend.text = element_text(angle = 90),
        panel.grid.major = element_line(color = "#b4b4b4",linewidth = .1),
        # panel.grid.minor = element_line(color = "black", size = .5),
        # panel.border = element_blank(),
        axis.line = element_line(color = "black")
  )+
  guides(fill = F)+  #隐藏代表物种颜色的图例
  scale_y_discrete(limits = unique(data_3$genefamily),position = "left")+
  scale_x_discrete(position = "bottom");p3
########################################################################################################################
########################################################################################################################
##加权柱状图
data_4 <- read_excel("~/R  project/207spe_ctdg/gene_family_counts.xlsx",sheet="Sheet2")
taxon_order <- rev(c("Angiosperms","Gymnosperms","Ferns","Lycophytes","Bryophytes","Zygnematophyceae",
                     "Streptophyte algae","Chlorophyta","Prasinodermophyta"))
# 确保 taxon 列是字符型
data_4$taxon <- as.character(data_4$taxon)

# 计算每个基因家族的权重
gene_families <- setdiff(names(data_4), "taxon")
weights <- 1 / sapply(data_4[gene_families], max, na.rm = TRUE)

# 计算加权平均
weighted_counts <- data_4 %>%
  mutate(across(all_of(gene_families), ~ .x * weights[cur_column()])) %>%
  rowwise() %>%
  mutate(avg_weighted_count = mean(c_across(all_of(gene_families)), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(avg_weighted_count = avg_weighted_count * mean(sapply(data_4[gene_families], max, na.rm = TRUE)))

# 计算每个类群的加权平均
taxon_avg <- weighted_counts %>%
  group_by(taxon) %>%
  summarise(avg_count = mean(avg_weighted_count, na.rm = TRUE),
            species_count = n())
# ##对taxon_avg进行排序，以便后面添加柱形图的边界信息
# taxon_avg$taxon <- factor(taxon_avg$taxon, levels = taxon_order)
# taxon_avg <- taxon_avg[order(taxon_avg$taxon), ]
# ##添加样本数量作为柱形图的宽度
# taxon_avg$right <- cumsum(taxon_avg$species_count) 
# taxon_avg$left <- taxon_avg$right - taxon_avg$species_count
# ##将上面得到的柱形图边界数据合并到点图数据中
# weighted_counts <- merge(weighted_counts, taxon_avg[, c("taxon", "left", "right")], by = "taxon", all.x = TRUE)



p4 <- ggplot() +
  # geom_rect(data = taxon_avg,aes(xmin = left, xmax = right, ymax = avg_count,ymin=0,fill = taxon),
  #           alpha=0.6)+ ## 柱子宽度反映物种数量
  geom_bar(data = taxon_avg,
           aes(x = factor(taxon, levels = taxon_order), y = avg_count, fill = taxon),
           stat = "identity", alpha = 0.7) +
  geom_point(data = weighted_counts,
             aes(x = factor(taxon, levels = taxon_order), y = avg_weighted_count),
             color = "black", alpha = 0.5, size = .1,
             position = position_jitter(width = 0.6, height = 0)) +
  geom_text(data = taxon_avg, 
            aes(x = factor(taxon, levels = taxon_order),  
                y = avg_count, 
                label = round(avg_count, 2)),  # 显示每个柱子的值，保留两位小数
            vjust = -1, hjust = 0.5, size = 3, color = "black") +  # 调整文本位置
  # geom_text(data = taxon_avg, 
  #           aes(x = (left + right) / 2, y = 0, label = paste("n =", species_count)), 
  #           vjust = 1.5, hjust = 0.5, size = 3, color = "black")+ ## 柱子宽度反映物种数量
  geom_text(data = taxon_avg, 
            aes(x = factor(taxon, levels = taxon_order),  y = 0, label = paste("n =", species_count)), 
            vjust = 1.5, hjust = 0.5, size = 4, color = "black")+
  scale_fill_manual(values = colors_2) +
  guides(fill=F)+
  theme_minimal() +
  theme(
        axis.text.x =element_blank(),
        # axis.text.x = element_text(angle = 30,hjust=1,size = 9,color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8,color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black",size = .3),
        axis.line.x = element_line(color = "black",size = .3),
        axis.ticks.x = element_line(color = "black",size = .3),
        axis.ticks.y = element_line(color = "black",size = .3),
        legend.position = "right") +
  labs(y = "Weighted average gene family size",
       title = "",
       fill = "Taxon Group") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA));p4

########################################################################################################################
########################################################################################################################
###家族扩增热图
data_5 <- read_excel("~/R  project/207spe_ctdg/gene_family_counts.xlsx",sheet="size change")
data_5$LatinName <- factor(data_5$LatinName,levels = rev(unique(data_5$LatinName)))
data_5 <- gather(data_5,genefamily,count,-LatinName)
data_5$color <- with(data_5, 
                    ifelse(count== "B", "#e6f5c9", 
                    ifelse(count == 0, "ns", 
                    ifelse(count > 1, "#fb9a99", 
                    ifelse(count > 0 & count <= 1, "#a6cee3", NA)))))

p5 <- ggplot(data_5, aes(x = LatinName, y = genefamily, fill = color)) +
  geom_tile(color = "#b4b4b4") +
  # scale_fill_identity() +
  scale_fill_manual(values = c("ns" = "white", "#e6f5c9" = "#e6f5c9", "#a6cee3" = "#a6cee3", "#fb9a99" = "#fb9a99"),
                    name = "expansion and contraction",
                    labels = c( "#e6f5c9" = "Birth", "#a6cee3" = "contraction", "#fb9a99" = "expansion")
                    # guide = guide_legend(title.theme = element_text(angle = 90))
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size =8,colour = "black"),
    axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.title.x = element_blank(),#不显示x轴名称
    axis.title.y = element_blank(),
    legend.position = "right",
    scale_y_discrete(limits = unique(data_4$LatinName),position = "left"),
    scale_x_discrete(position = "top"),
    #axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(color = "black", size = 0.5),
    legend.text = element_text(angle = 0));p5

########################################################################################################################
########################################################################################################################
# p4 <- plot_grid(p1, p2,p3, ncol = 3,rel_heights = c(1,1,1),rel_widths = c(1,1,1))
p <- plot_grid(p4,p3,nrow = 2,rel_heights = c(5,16));p
# p <- plot_grid(p4,p3,p5,nrow = 3,rel_heights = c(5,16,5));p
setwd("~/Desktop/project/CTDG_of_207_spe/draft_fig/")
# ggsave("gene_family_8.26*9.pdf", plot = p, width = 8.26, height = 9, units = "in")
ggsave("gene_family_size_change_8.26*3.2.pdf", plot = p5, width = 8.26, height = 3.2, units = "in")
