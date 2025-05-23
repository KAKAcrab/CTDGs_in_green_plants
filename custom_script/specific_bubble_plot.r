#!/usr/bin/env Rscript

# 2OG-FeII_Oxy家族基因分布气泡图绘制脚本
# 作者: [您的姓名]
# 日期: 2025年5月
# 用途: 绘制2OG-FeII_Oxy家族基因在220个植物物种中的分布情况

# 加载所需的R包 / Load required packages
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(optparse)

# 定义命令行参数 / Define command line arguments
option_list <- list(
  make_option(c("-s", "--species_list"), type="character", default=NULL, 
              help="物种列表文件路径（220_spe_list_low2high_latin_abbr.csv）", metavar="character"),
  make_option(c("-t", "--total_genes"), type="character", default=NULL, 
              help="总基因ID文件路径（220_total_gene_ID.csv）", metavar="character"),
  make_option(c("-c", "--clustered_genes"), type="character", default=NULL, 
              help="成簇基因ID文件路径（220_clustered_gene_ID.csv）", metavar="character"),
  make_option(c("-g", "--gene_clusters"), type="character", default=NULL, 
              help="基因簇ID文件路径（220_CTDGs_ID.csv）", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="output", 
              help="输出目录", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 验证输入参数 / Validate input parameters
if (is.null(opt$species_list) || is.null(opt$total_genes) || 
    is.null(opt$clustered_genes) || is.null(opt$gene_clusters)) {
  stop("请提供所有必需的输入文件路径！\nPlease provide all required input file paths!")
}

# 创建输出目录 / Create output directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# 定义三种指标的颜色 / Define colors for three metrics
metric_colors <- c(
  "total_genes" = "#fbb4ae",        # 蓝色 - 总基因数 / Blue - Total genes
  "clustered_genes" = "#b3cde3",    # 紫色 - 成簇基因数 / Purple - Clustered genes
  "gene_clusters" = "#ccebc5"       # 浅紫色 - 基因簇数 / Light purple - Gene clusters
)

# 定义指标的英文标签 / Define English labels for metrics
metric_labels <- c(
  "total_genes" = "Total genes number",
  "clustered_genes" = "Clustered genes number", 
  "gene_clusters" = "CTDGs number"
)

# 1. 读取物种列表 / Read species list
cat("读取物种列表文件...\nReading species list file...\n")
species_list <- read.csv(opt$species_list, header = TRUE, stringsAsFactors = FALSE)
# 假设第一列是拉丁名，第二列是缩写
# Assume first column is Latin name, second column is abbreviation
colnames(species_list)[1:2] <- c("latin_name", "species_abbr")

cat("共找到", nrow(species_list), "个物种\nFound", nrow(species_list), "species\n")

# 2. 读取基因数据 / Read gene data
cat("读取基因数据文件...\nReading gene data files...\n")

# 读取总基因数据 / Read total gene data
total_genes_data <- read.csv(opt$total_genes, header = TRUE, stringsAsFactors = FALSE)
colnames(total_genes_data) <- c("species_abbr", "gene_id")

# 读取成簇基因数据 / Read clustered gene data  
clustered_genes_data <- read.csv(opt$clustered_genes, header = TRUE, stringsAsFactors = FALSE)
colnames(clustered_genes_data) <- c("species_abbr", "gene_id")

# 读取基因簇数据 / Read gene cluster data
gene_clusters_data <- read.csv(opt$gene_clusters, header = TRUE, stringsAsFactors = FALSE)
colnames(gene_clusters_data) <- c("species_abbr", "cluster_id")

# 3. 统计每个物种的基因数量 / Count genes for each species
cat("统计每个物种的基因数量...\nCounting genes for each species...\n")

# 统计总基因数 / Count total genes
total_gene_counts <- total_genes_data %>%
  group_by(species_abbr) %>%
  summarise(total_genes = n(), .groups = 'drop')

# 统计成簇基因数 / Count clustered genes
clustered_gene_counts <- clustered_genes_data %>%
  group_by(species_abbr) %>%
  summarise(clustered_genes = n(), .groups = 'drop')

# 统计基因簇数 / Count gene clusters
gene_cluster_counts <- gene_clusters_data %>%
  group_by(species_abbr) %>%
  summarise(gene_clusters = n(), .groups = 'drop')

# 4. 合并数据并补充缺失物种 / Merge data and add missing species
cat("合并数据并补充缺失物种的零值...\nMerging data and adding zero values for missing species...\n")

# 创建完整的物种数据框 / Create complete species data frame
complete_data <- species_list %>%
  left_join(total_gene_counts, by = "species_abbr") %>%
  left_join(clustered_gene_counts, by = "species_abbr") %>%
  left_join(gene_cluster_counts, by = "species_abbr") %>%
  # 将NA替换为0 / Replace NA with 0
  mutate(
    total_genes = ifelse(is.na(total_genes), 0, total_genes),
    clustered_genes = ifelse(is.na(clustered_genes), 0, clustered_genes),
    gene_clusters = ifelse(is.na(gene_clusters), 0, gene_clusters)
  )

# 5. 将数据转换为长格式用于绘图 / Convert data to long format for plotting
cat("准备绘图数据...\nPreparing plotting data...\n")

plot_data <- complete_data %>%
  select(latin_name, species_abbr, total_genes, clustered_genes, gene_clusters) %>%
  pivot_longer(
    cols = c(total_genes, clustered_genes, gene_clusters),
    names_to = "metric",
    values_to = "count"
  ) %>%
  # 设置因子水平以控制显示顺序 / Set factor levels to control display order
  mutate(
    latin_name = factor(latin_name, levels = species_list$latin_name),
    metric = factor(metric, levels = c("total_genes", "clustered_genes", "gene_clusters"))
  )

# 添加指标标签 / Add metric labels
plot_data$metric_label <- metric_labels[as.character(plot_data$metric)]
plot_data$metric_label <- factor(plot_data$metric_label, 
                                 levels = c(metric_labels[["total_genes"]], 
                                            metric_labels[["clustered_genes"]], 
                                            metric_labels[["gene_clusters"]]))

# 6. 创建气泡图 / Create bubble plot
cat("绘制气泡图...\nCreating bubble plot...\n")

# 计算合适的气泡大小范围 / Calculate appropriate bubble size range
max_count <- max(plot_data$count)
cat("最大基因数量:", max_count, "\nMaximum gene count:", max_count, "\n")

# 创建气泡图 / Create bubble plot
bubble_plot <- ggplot(plot_data, aes(x = metric_label, y = latin_name)) +
  geom_point(aes(size = count, fill = metric), 
             shape = 21, color = "black", alpha = 0.8, stroke = 0.3) +
  # 设置气泡大小映射（不填充颜色）/ Set bubble size mapping (no fill color)
  scale_size_area(max_size = 8, 
                  name = "Gene count",
                  breaks = c(0, 5, 10, 20, 50, 100),
                  labels = c("0", "5", "10", "20", "50", "100+"),
                  guide = guide_legend(override.aes = list(color = "black", 
                                                           fill = "white", 
                                                           stroke = 0.5))) +
  # 设置颜色映射（不显示图例）/ Set color mapping (no legend)
  scale_fill_manual(values = metric_colors, guide = "none") +
  # 设置坐标轴标签 / Set axis labels
  labs(x = "", y = "") +
  # 设置主题 / Set theme
  theme_minimal() +
  theme(
    # 设置坐标轴文字 / Set axis text
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, 
                               size = 8, color = "black", face = "bold"),
    axis.text.y = element_text(size = 4, color = "black", face = "italic"),
    # 移除网格线 / Remove grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 添加面板边框 / Add panel border
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # 设置图例位置和样式 / Set legend position and style
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(l = 5),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.title.position = "top",
    # 设置图边距 / Set plot margins
    plot.margin = margin(t = 0, r = 3, b = 3, l = 50)
  ) +
  # 设置坐标轴位置 / Set axis position
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(plot_data$latin_name)))

# 7. 保存图表 / Save plot
cat("保存图表...\nSaving plot...\n")

# 构建完整的输出路径 / Build complete output path
output_pdf <- file.path(opt$output_dir, "2OG-FeII_Oxy_gene_distribution_bubble_plot.pdf")

# 保存为PDF格式 / Save as PDF
ggsave(filename = output_pdf, 
       plot = bubble_plot, 
       width = 8.3, height = 11.7, units = "in")

# 验证文件是否成功保存 / Verify file was saved successfully
if (file.exists(output_pdf)) {
  cat("PDF文件已成功保存到:", output_pdf, "\nPDF file successfully saved to:", output_pdf, "\n")
} else {
  cat("警告：PDF文件保存失败！\nWarning: PDF file save failed!\n")
}

# 8. 保存统计数据 / Save statistics
cat("保存统计数据...\nSaving statistics...\n")

# 保存完整数据 / Save complete data
write.csv(complete_data, 
          file = file.path(opt$output_dir, "2OG-FeII_Oxy_gene_statistics.csv"), 
          row.names = FALSE)

# 保存汇总统计 / Save summary statistics
summary_stats <- complete_data %>%
  summarise(
    total_species = n(),
    species_with_total_genes = sum(total_genes > 0),
    species_with_clustered_genes = sum(clustered_genes > 0),
    species_with_gene_clusters = sum(gene_clusters > 0),
    max_total_genes = max(total_genes),
    max_clustered_genes = max(clustered_genes),
    max_gene_clusters = max(gene_clusters),
    mean_total_genes = round(mean(total_genes), 2),
    mean_clustered_genes = round(mean(clustered_genes), 2),
    mean_gene_clusters = round(mean(gene_clusters), 2)
  )

write.csv(summary_stats, 
          file = file.path(opt$output_dir, "2OG-FeII_Oxy_gene_summary_statistics.csv"), 
          row.names = FALSE)

# 9. 输出运行信息 / Output run information
cat("\n=== 运行完成 / Run Completed ===\n")
cat("处理的物种数量 / Number of species processed:", nrow(complete_data), "\n")
cat("输出文件保存在 / Output files saved in:", opt$output_dir, "\n")
cat("- 气泡图PDF / Bubble plot PDF: 2OG-FeII_Oxy_gene_distribution_bubble_plot.pdf\n")
cat("- 统计数据 / Statistics: 2OG-FeII_Oxy_gene_statistics.csv\n")
cat("- 汇总统计 / Summary statistics: 2OG-FeII_Oxy_gene_summary_statistics.csv\n")

# 显示部分数据预览 / Show data preview
cat("\n=== 数据预览 / Data Preview ===\n")
print(head(complete_data, 10))

cat("\n=== 汇总统计 / Summary Statistics ===\n")
print(summary_stats)

cat("\n脚本执行完成！\nScript execution completed!\n")