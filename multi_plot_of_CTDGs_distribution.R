library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)
library(viridisLite)
library(viridis)
library(hrbrthemes)
library(ggpubr)
library(gridExtra)
library(patchwork)
library(writexl)
##
setwd("~/Desktop/project/CTDG_of_207_spe/draft_fig/")
light_colors <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462",
                  "#b3de69","#a6cee3","#1f78b4","#fb9a99","#6a3d9a","#b15928","#e31a1c")

multi_data <- read.csv("~/R  project/207spe_ctdg/169_distribution_of_CTDGs_in_chromosome.csv",
                                     header = TRUE,sep =',')
multi_data <- multi_data %>%
  mutate(
    chr_len = gsub("[^0-9]", "", chr_len),  # 去除非数字字符
    CTDG_num = gsub("[^0-9]", "", CTDG_num) # 去除非数字字符
  )  
multi_data <- multi_data %>%
  mutate(
    chr_len = as.numeric(chr_len),
    CTDG_num = as.numeric(CTDG_num)
  )
multi_data <- multi_data[,c("species","chr","chr_num","chr_len","CTDG_num")]
multi_data$chr_len_mb <- multi_data$chr_len / 1e6
multi_data$species <- factor(multi_data$species, levels = unique(multi_data$species))
head(multi_data)
# 计算每个物种的总CTDG数量和平均染色体长度
species_summary <- multi_data %>%
  group_by(species) %>%
  summarise(chr_count = dplyr::first(chr_num))
# 计算每个物种的染色体长度和CTDG数量的相关系数
# 检查标准差是否为零，并且在计算相关性之前过滤掉零标准差的组
cor_data <- multi_data %>%
  group_by(species) %>%
  summarise(
    correlation = ifelse(sd(chr_len) > 0 & sd(CTDG_num) > 0, cor(chr_len, CTDG_num,method = "pearson"), NA)
  ) %>%
  mutate(color = ifelse(correlation > 0, "#fb8072", "#a6cee3"))

# 去掉 NA 值
cor_data <- cor_data %>%
  filter(!is.na(correlation))

# # 将 cor_data 合并到 multi_data，依据 species 列进行合并
# multi_data_with_correlation <- multi_data %>%
#   left_join(cor_data, by = "species")
# write_xlsx(multi_data_with_correlation, "CTDG_distribution_with_correlation.xlsx")

a <- cor_data %>% filter(correlation < 0) %>% nrow();a


# 绘制条形图 - 第一列子图
p1 <- ggplot(species_summary, aes(x = chr_count,y = species)) +
  geom_bar(aes(x = chr_count,y = species),stat = "identity", position = "identity", fill="#80b1d3") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "black",size = .2,linetype = "dashed"),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = .5),
        axis.text.x=element_text(angle = 0,hjust = 0.5,size = 8,color = "black"),
        axis.title.x = element_text(color = "black",size = 10),
        # axis.text.y = element_blank(),
        # axis.text.x=element_blank(),
        axis.text.y = element_text(size = 4,color = "black"),
        # axis.ticks.y = element_blank(),
        axis.ticks.x= element_line(color = "black"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(0, 0, 0, 0),
        # strip.text = element_text(size = 8),#调整每个分组标题的字体大小
        legend.position = "top"
  )+
  scale_y_discrete(limits = rev(levels(species_summary$species)))+
  labs(title = "", x = "Number of chromosomes", y = "");p1

# 绘制染色体长度的箱线图 - 第二列子图
p2 <- ggplot(multi_data, aes(y = species, x = chr_len_mb)) +
  geom_boxplot(cex = .1,
               outlier.size =.5,
               show.legend = T,
               color="#e31a1c") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "black",size = .2,linetype = "dashed"),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = .5),
        axis.text.x=element_text(angle = 0,hjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(color = "black",size = 10),
        axis.ticks.x= element_line(color = "black"),
        # axis.text.y = element_text(size = 5,color = "black"),
        plot.title = element_text(size = 8, face = "bold", hjust = .5),
        strip.text = element_text(size = 8),#调整每个分组标题的字体大小
        legend.position = "top"
  )+
  labs(title = "", x = "Chromosome length(Mb)", y = "");p2

# 绘制CTDG数量的箱线图 - 第三列子图
p3 <- ggplot(multi_data, aes(y = species, x = CTDG_num)) +
  geom_boxplot(cex = .1,
               outlier.size =.5,
               show.legend = T,
               color="#6a3d9a") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "black",size = .2,linetype = "dashed"),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = .5),
        axis.text.x=element_text(angle = 0,hjust = 0.5,size = 8,color = "black"),
        axis.title.x = element_text(color = "black",size = 10),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size = 5,color = "black"),
        axis.ticks.x= element_line(color = "black"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        strip.text = element_text(size = 8),#调整每个分组标题的字体大小
        legend.position = "top"
  )+
  labs(title = "", x = "Number of CTDGs", y = "");p3

# 绘制Pearson相关系数的条形图 - 第四列子图
p4 <- ggplot(cor_data, aes(y = species, x = correlation, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "black",size = .2,linetype = "dashed"),
        axis.line.x=element_blank(),
        axis.line.y=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = .5),
        axis.text.x=element_text(angle = 0,hjust = 0.5,size = 8,color = "black"),
        axis.title.x = element_text(color = "black",size = 10,hjust = .5),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size = 5,color = "black"),
        axis.ticks.x= element_line(color = "black"),
        axis.title.y.right = element_text(size = 8,color = "black"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        strip.text = element_text(size = 8),#调整每个分组标题的字体大小
        legend.position = "top"
  )+
  labs(title = "", x = "Pearson correlation\ncoefficient", y = "");p4

# p0 <- ggplot(data.frame(species = levels(multi_data$species), y = 1:length(levels(multi_data$species))),
#                         aes(x = 0, y = y, label = species)) +
#       # geom_text(hjust = 1, size = 3) +
#   scale_y_reverse(breaks = seq_along(levels(multi_data$species)), labels = levels(multi_data$species)) +
#       theme_void() +
#   coord_cartesian(clip = "off")+
#       theme(
#         plot.margin = margin(5.5, 0, 5.5, 5.5),
#         # axis.text.y = element_blank()
#         axis.text.y = element_text(size = 6,color = "black",hjust = 1)
#         );p0
p <- ( p1 + p2 + p3 + p4) + 
  plot_layout(ncol = 4, widths = c( 1, 1, 1, 1),guides = "collect");p
ggsave("Figure 5.pdf",
       plot = p, width = 8.26, height = 11.69, units = "in")

