#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(optparse)

# 设置命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="输入CSV文件路径", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="output", 
              help="输出目录 [默认= %default]", metavar="character"),
  make_option(c("-t", "--top_n"), type="integer", default=30, 
              help="展示变化最大的前N个基因 [默认= %default]", metavar="integer"),
  make_option(c("-m", "--min_species"), type="integer", default=10, 
              help="基因至少在多少个物种中有数据才会被考虑 [默认= %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 验证输入参数
if (is.null(opt$input)) {
  stop("输入CSV文件路径 (-i 或 --input)")
}

# 创建输出目录
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# 读取数据
cat("正在读取数据文件:", opt$input, "\n")
data <- read.csv(opt$input, stringsAsFactors = FALSE, check.names = FALSE)

# 确认数据包含物种列
if (!"Species" %in% colnames(data)) {
  stop("CSV文件必须包含'Species'列")
}

# 从列名中提取基因信息
cat("正在解析基因列...\n")
# 获取成簇基因列
gene_columns <- grep("list_[1-7]_AT.*[^_clustered]$", colnames(data), value = TRUE)
# 获取完整数量基因列
clustered_columns <- grep("list_[1-7]_AT.*_clustered$", colnames(data), value = TRUE)

# 检查是否找到基因列
if (length(gene_columns) == 0 || length(clustered_columns) == 0) {
  stop("无法从CSV中找到基因和成簇基因列")
}

cat("找到", length(gene_columns), "个基因列和", length(clustered_columns), "个成簇基因列\n")

# 检查基因与成簇基因列的匹配
base_gene_names <- sub("_clustered$", "", clustered_columns)
if (!all(base_gene_names %in% gene_columns)) {
  warning("部分成簇基因列没有对应的基因总数列")
}

# 创建比例数据框
cat("计算成簇比例...\n")
ratio_data <- data.frame(Species = data$Species)

# 计算每个基因的成簇比例
for (gene in base_gene_names) {
  if (gene %in% gene_columns) {
    clustered_col <- paste0(gene, "_clustered")
    
    # 计算比例，处理分母为0的情况
    ratio_data[[gene]] <- ifelse(data[[gene]] > 0, 
                                 data[[clustered_col]] / data[[gene]], 
                                 NA)
  }
}

# 删除缺失值过多的基因
valid_columns <- c("Species")
for (col in colnames(ratio_data)[-1]) {
  valid_count <- sum(!is.na(ratio_data[[col]]))
  if (valid_count >= opt$min_species) {
    valid_columns <- c(valid_columns, col)
  } else {
    cat("移除基因", col, "（仅在", valid_count, "个物种中有数据）\n")
  }
}
ratio_data <- ratio_data[, valid_columns, drop = FALSE]

cat("成簇比例计算完成，保留了", length(valid_columns)-1, "个有效基因\n")

# 计算每个基因的变化度量
cat("分析基因变化趋势...\n")
gene_variation <- data.frame(gene = character(), 
                             variation = numeric(), # 变化幅度，最大值与最小值的差，反映极端变化
                             cv = numeric(),  # 变异系数，标准差除以均值，提供相对变异度量
                             mean_ratio = numeric(), # 平均比例，反映基因的总体成簇倾向
                             stringsAsFactors = FALSE)

for (gene in colnames(ratio_data)[-1]) {
  values <- ratio_data[[gene]]
  values <- values[!is.na(values)]
  
  if (length(values) >= opt$min_species) {
    variation <- max(values) - min(values)
    cv <- sd(values) / mean(values)  # 变异系数
    mean_val <- mean(values)
    
    gene_variation <- rbind(gene_variation, 
                            data.frame(gene = gene, 
                                       variation = variation,
                                       cv = cv,
                                       mean_ratio = mean_val,
                                       stringsAsFactors = FALSE))
  }
}

# 按变化幅度排序
gene_variation <- gene_variation[order(-gene_variation$variation), ]

# 提取前N个变化最大的基因
top_n_genes <- min(opt$top_n, nrow(gene_variation))
top_varying_genes <- head(gene_variation, top_n_genes)
cat("已识别出变化最大的", nrow(top_varying_genes), "个基因\n")

# 准备热图数据
heatmap_data <- ratio_data
rownames(heatmap_data) <- heatmap_data$Species
heatmap_data$Species <- NULL

# 预处理热图数据 - 根据变化率对基因排序
ordered_genes <- gene_variation$gene
heatmap_data_ordered <- heatmap_data[, ordered_genes, drop = FALSE]

# 准备top基因热图数据
top_genes <- top_varying_genes$gene
heatmap_data_top <- heatmap_data[, top_genes, drop = FALSE]

# 保存变化最大的基因列表
write.csv(gene_variation, file.path(opt$output_dir, "gene_variation_ranking.csv"), row.names = FALSE)
cat("变化最大的基因排名已保存到:", file.path(opt$output_dir, "gene_variation_ranking.csv"), "\n")

# 提取性状和基因ID信息
extract_info <- function(gene_name) {
  parts <- strsplit(gene_name, "_")[[1]]
  trait <- paste0(parts[1], "_", parts[2])
  gene_id <- parts[3]
  return(c(trait = trait, gene_id = gene_id))
}

gene_info <- t(sapply(colnames(heatmap_data_top), extract_info))
gene_labels <- paste0(gene_info[, "trait"], "_", gene_info[, "gene_id"])

# 绘制热图
cat("生成热图...\n")

# 设置颜色渐变
my_colors <- colorRampPalette(c("white", "yellow", "red"))(100)

# 所有基因的热图
pdf(file.path(opt$output_dir, "cluster_ratio_heatmap_all_genes.pdf"), width = 20, height = 20)
pheatmap(heatmap_data_ordered,
         color = my_colors,
         cluster_rows = FALSE,  # 保持物种原有顺序
         cluster_cols = TRUE,   # 基因聚类
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Clustered proportional heat map of all genes",
         fontsize = 7,
         fontsize_col = 6)
dev.off()

# Top变化基因的热图
pdf(file.path(opt$output_dir, "cluster_ratio_heatmap_top_genes.pdf"), width = 14, height = 20)
pheatmap(heatmap_data_top,
         color = my_colors,
         cluster_rows = FALSE,  # 保持物种原有顺序
         cluster_cols = TRUE,   # 基因聚类
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = paste0("Top ", nrow(top_varying_genes), "clustered propotional genes heat map"),
         fontsize = 7,
         fontsize_col = 8)
dev.off()

# 为折线图准备数据 - 前10个变化最大的基因
top10_genes <- head(top_varying_genes, 10)$gene
line_data <- melt(ratio_data[, c("Species", top10_genes)], 
                  id.vars = "Species", 
                  variable.name = "Gene", 
                  value.name = "Ratio")

# 提取基因的trait信息
line_data$Trait <- sapply(as.character(line_data$Gene), function(g) {
  paste0(strsplit(g, "_")[[1]][1:2], collapse = "_")
})

# 绘制折线图
cat("生成折线图...\n")
p <- ggplot(line_data, aes(x = Species, y = Ratio, color = Gene, group = Gene)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Trait, scales = "free_y",nrow = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)))+
  scale_x_discrete(limits = rev(unique(line_data$Species)))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "top",
        legend.direction = "horizontal") +
  labs(title = "Cluster ratio trends of top 10 genes",
       x = "Species(aquatic --> terrestrial)",
       y = "Cluster ratio",
       color = "Gene")

# 保存折线图
ggsave(file.path(opt$output_dir, "cluster_ratio_line_plot.pdf"), p, width = 14, height = 10)

# 计算进化趋势
cat("分析进化趋势...\n")
# 创建物种进化水平的数值表示（按照CSV中的顺序，假定是从低等到高等）
species_evolution_rank <- seq_len(nrow(ratio_data))
names(species_evolution_rank) <- ratio_data$Species

# 计算每个基因的进化趋势（使用相关系数）
evolution_trend <- data.frame(gene = character(), 
                              correlation = numeric(), 
                              p_value = numeric(), 
                              stringsAsFactors = FALSE)

for (gene in colnames(ratio_data)[-1]) {
  gene_data <- ratio_data[[gene]]
  valid_idx <- !is.na(gene_data)
  
  if (sum(valid_idx) >= opt$min_species) {  # 至少需要一定数量的有效数据点
    cor_test <- cor.test(species_evolution_rank[valid_idx], 
                         gene_data[valid_idx],
                         method = "spearman")
    
    evolution_trend <- rbind(evolution_trend, 
                             data.frame(gene = gene,
                                        correlation = cor_test$estimate,
                                        p_value = cor_test$p.value,
                                        stringsAsFactors = FALSE))
  }
}

# 添加多重测试校正
evolution_trend$p_adj <- p.adjust(evolution_trend$p_value, method = "BH")

# 排序并找出显著相关的基因
evolution_trend <- evolution_trend[order(evolution_trend$p_adj), ]
significant_genes <- evolution_trend[evolution_trend$p_adj < 0.05, ]

# 保存进化趋势分析结果
write.csv(evolution_trend, file.path(opt$output_dir, "evolution_trend_analysis.csv"), row.names = FALSE)
cat("进化趋势分析结果已保存到:", file.path(opt$output_dir, "evolution_trend_analysis.csv"), "\n")

# 可视化显著的进化趋势
if (nrow(significant_genes) > 0) {
  # 为折线图准备数据
  sig_genes <- head(significant_genes$gene, 20)  # 最多取前20个显著基因
  trend_data <- melt(ratio_data[, c("Species", sig_genes)], 
                     id.vars = "Species", 
                     variable.name = "Gene", 
                     value.name = "Ratio")
  
  # 提取基因的trait信息
  trend_data$Trait <- sapply(as.character(trend_data$Gene), function(g) {
    paste0(strsplit(g, "_")[[1]][1:2], collapse = "_")
  })
  
  # 根据相关系数为基因增加方向信息
  gene_direction <- setNames(significant_genes$correlation, significant_genes$gene)
  trend_data$Direction <- ifelse(gene_direction[as.character(trend_data$Gene)] > 0, 
                                 "增加趋势", "减少趋势")
  
  # 绘制折线图，按进化趋势分组
  p <- ggplot(trend_data, aes(x = Species, y = Ratio, color = Direction, group = Gene)) +
    geom_line(alpha = 0.7) +
    scale_color_manual(values = c("Increasing trend" = "red", "Decreasing trend" = "blue")) +
    facet_wrap(~Trait, scales = "free_y",nrow = 6) +
    scale_x_discrete(limits = rev(unique(line_data$Species)))+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
    labs(title = "Significant evolution trends of cluster ratio",
         subtitle = paste0("Significant genes number: ", nrow(significant_genes)),
         x = "Species(aquatic --> terrestrial)",
         y = "Cluster ratio",
         color = "Trends")
  
  # 保存进化趋势折线图
  ggsave(file.path(opt$output_dir, "evolution_trend_line_plot.pdf"), p, width = 14, height = 20)
  cat("进化趋势图已保存到:", file.path(opt$output_dir, "evolution_trend_line_plot.pdf"), "\n")
  
  # 创建候选基因列表（变化大且与进化相关的基因）
  candidate_genes <- merge(top_varying_genes, significant_genes, by = "gene")
  candidate_genes <- candidate_genes[order(-abs(candidate_genes$correlation)), ]
  write.csv(candidate_genes, file.path(opt$output_dir, "candidate_genes_for_phylogeny.csv"), row.names = FALSE)
  cat("最佳系统发育分析候选基因已保存到:", file.path(opt$output_dir, "candidate_genes_for_phylogeny.csv"), "\n")
}

cat("\n分析完成! 所有结果已保存到目录:", opt$output_dir, "\n")