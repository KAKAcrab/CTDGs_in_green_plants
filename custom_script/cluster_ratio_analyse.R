#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(optparse)
library(writexl)  # 用于写入xlsx文件

# 设置命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="输入CSV文件路径", metavar="character"),
  make_option(c("-g", "--group_map"), type="character", default=NULL, 
              help="类群物种映射文件路径 (group2spe.csv)", metavar="character"),
  make_option(c("-n", "--name_map"), type="character", default=NULL, 
              help="基因ID到家族名称映射文件路径 (gene_family_name2_gene_ID.csv)", metavar="character"),
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
if (is.null(opt$input) || is.null(opt$group_map) || is.null(opt$name_map)) {
  stop("请提供所有必需的输入文件路径: 成簇率数据(-i), 类群映射(-g), 基因名称映射(-n)")
}

# 创建输出目录
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# 读取数据
cat("正在读取数据文件:", opt$input, "\n")
data <- read.csv(opt$input, stringsAsFactors = FALSE, check.names = FALSE)

# 读取类群物种映射
cat("正在读取类群物种映射:", opt$group_map, "\n")
group2spe <- read.csv(opt$group_map, header = FALSE, 
                      col.names = c("Group", "Species"),
                      stringsAsFactors = FALSE)

# 读取基因家族映射
cat("正在读取基因名称映射:", opt$name_map, "\n")
gene_family_map <- read.csv(opt$name_map, header = FALSE, 
                            col.names = c("Family", "GeneID"),
                            stringsAsFactors = FALSE)

# 创建基因ID到家族名称的映射字典
gene_id_to_family <- setNames(gene_family_map$Family, gene_family_map$GeneID)

cat("读取到", nrow(gene_family_map), "个基因家族映射\n")
if(nrow(gene_family_map) > 0) {
  cat("示例映射：", gene_family_map$GeneID[1], "->", gene_family_map$Family[1], "\n")
}

# 确认数据包含物种列
if (!"Species" %in% colnames(data)) {
  stop("CSV文件必须包含'Species'列")
}

# 从列名中提取基因信息
cat("正在解析基因列...\n")
# 获取成簇基因列
gene_columns <- grep("^list_[1-7]_.*[^_clustered]$", colnames(data), value = TRUE)
# 获取完整数量基因列
clustered_columns <- grep("^list_[1-7]_.*_clustered$", colnames(data), value = TRUE)

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

# 功能：将基因ID转换为家族名称
convert_gene_names <- function(names_vector) {
  converted_names <- names_vector
  
  for (i in seq_along(names_vector)) {
    gene_name <- names_vector[i]
    # 跳过Species列
    if (gene_name == "Species" || gene_name == "Group") {
      next
    }
    
    # 检查完整的基因名是否在映射中
    if (gene_name %in% names(gene_id_to_family)) {
      converted_names[i] <- gene_id_to_family[gene_name]
      cat("替换基因名:", gene_name, "->", converted_names[i], "\n")
    } else {
      # 尝试提取基因ID部分
      parts <- strsplit(gene_name, "_")[[1]]
      if (length(parts) >= 3) {
        # 重构完整的基因ID
        trait_prefix <- paste0(parts[1], "_", parts[2])
        gene_id <- parts[3]
        full_id <- paste0(trait_prefix, "_", gene_id)
        
        if (full_id %in% names(gene_id_to_family)) {
          converted_names[i] <- gene_id_to_family[full_id]
          cat("替换基因名(部分匹配):", gene_name, "->", converted_names[i], "\n")
        }
      }
    }
  }
  
  return(converted_names)
}

# 计算每个基因的变化度量
cat("分析基因变化趋势...\n")
gene_variation <- data.frame(gene = character(), 
                             variation = numeric(), 
                             cv = numeric(),  
                             mean_ratio = numeric(), 
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

# 转换基因ID为家族名称
gene_variation$family_name <- convert_gene_names(gene_variation$gene)

# 提取前N个变化最大的基因
top_n_genes <- min(opt$top_n, nrow(gene_variation))
top_varying_genes <- head(gene_variation, top_n_genes)
cat("已识别出变化最大的", nrow(top_varying_genes), "个基因\n")

# 保存变化最大的基因列表（使用家族名称）
write.csv(gene_variation[, c("family_name", "variation", "cv", "mean_ratio")], 
          file.path(opt$output_dir, "gene_variation_ranking.csv"), row.names = FALSE)
cat("变化最大的基因排名已保存到:", file.path(opt$output_dir, "gene_variation_ranking.csv"), "\n")

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

# 将基因ID转换为家族名称用于热图标签
colnames(heatmap_data_ordered) <- convert_gene_names(colnames(heatmap_data_ordered))
colnames(heatmap_data_top) <- convert_gene_names(colnames(heatmap_data_top))

# 绘制热图
cat("生成热图...\n")

# 设置颜色渐变
my_colors <- colorRampPalette(c("white", "yellow", "red"))(100)

# 所有基因的热图
heatmap_data_ordered[!is.finite(as.matrix(heatmap_data_ordered))] <- 0
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
heatmap_data_top[!is.finite(as.matrix(heatmap_data_top))] <- 0
pdf(file.path(opt$output_dir, "cluster_ratio_heatmap_top_genes.pdf"), width = 14, height = 20)
pheatmap(heatmap_data_top,
         color = my_colors,
         cluster_rows = FALSE,  # 保持物种原有顺序
         cluster_cols = TRUE,   # 基因聚类
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = paste0("Top ", nrow(top_varying_genes), " clustered propotional genes heat map"),
         fontsize = 7,
         fontsize_col = 8)
dev.off()

# 计算进化趋势
cat("分析进化趋势...\n")
# 创建物种进化水平的数值表示（按照CSV文件中的顺序，假定是从低等到高等）
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

# 转换基因ID为家族名称
evolution_trend$family_name <- convert_gene_names(evolution_trend$gene)

# 排序并找出显著相关的基因
evolution_trend <- evolution_trend[order(evolution_trend$p_adj), ]
significant_genes <- evolution_trend[evolution_trend$p_adj < 0.05, ]

# 保存进化趋势分析结果（使用家族名称）
write.csv(evolution_trend[, c("family_name", "correlation", "p_value", "p_adj")], 
          file.path(opt$output_dir, "evolution_trend_analysis.csv"), row.names = FALSE)
cat("进化趋势分析结果已保存到:", file.path(opt$output_dir, "evolution_trend_analysis.csv"), "\n")

# 创建候选基因列表（变化大且与进化相关的基因）
# 使用原始基因ID进行合并，然后在输出时使用家族名称
if (nrow(significant_genes) > 0) {
  candidate_genes <- merge(top_varying_genes, significant_genes, by = "gene")
  candidate_genes <- candidate_genes[order(-abs(candidate_genes$correlation)), ]
  
  # 输出使用家族名称的候选基因列表
  write.csv(candidate_genes[, c("family_name.x", "variation", "cv", "mean_ratio", 
                                "correlation", "p_value", "p_adj")], 
            file.path(opt$output_dir, "candidate_genes_for_phylogeny.csv"), 
            row.names = FALSE)
  
  # 重命名列以提高可读性
  colnames_mapping <- c(
    "family_name.x" = "family_name",
    "variation" = "variation",
    "cv" = "cv",
    "mean_ratio" = "mean_ratio",
    "correlation" = "correlation",
    "p_value" = "p_value",
    "p_adj" = "p_adj"
  )
  
  # 读取并重命名
  candidate_genes_csv <- read.csv(file.path(opt$output_dir, "candidate_genes_for_phylogeny.csv"))
  colnames(candidate_genes_csv) <- c("gene_family", "variation", "cv", "mean_ratio", 
                                     "correlation", "p_value", "p_adj")
  write.csv(candidate_genes_csv, 
            file.path(opt$output_dir, "candidate_genes_for_phylogeny.csv"), 
            row.names = FALSE)
  
  cat("最佳系统发育分析候选基因已保存到:", 
      file.path(opt$output_dir, "candidate_genes_for_phylogeny.csv"), "\n")
}

# 计算按类群聚合的成簇率（针对top_n个基因）
cat("计算类群平均成簇率...\n")

# 添加类群信息到ratio_data
ratio_data$Group <- sapply(ratio_data$Species, function(sp) {
  group_idx <- which(group2spe$Species == sp)
  if (length(group_idx) > 0) {
    return(group2spe$Group[group_idx[1]])
  } else {
    return(NA)
  }
})

# 移除没有对应类群的物种
ratio_data <- ratio_data[!is.na(ratio_data$Group), ]

# 类群顺序 - 从高等到低等植物
group_order <- c(
  "Angiosperms", 
  "Gymnosperms", 
  "Ferns",
  "Lycophytes", 
  "Bryophytes", 
  "Zygnematophyceae",  
  "Streptophyte_algae",  
  "Chlorophyta"
)

# 确保Group按指定顺序排序
ratio_data$Group <- factor(ratio_data$Group, levels = group_order)

# 计算每个类群每个基因的平均成簇率
group_ratios <- ratio_data %>%
  group_by(Group) %>%
  summarise(across(-Species, ~mean(.x, na.rm = TRUE)))

# 选择top_n个基因
top_genes_for_xlsx <- top_varying_genes$gene

# 子集出仅包含top_n个基因的数据
group_ratios_top <- group_ratios[, c("Group", top_genes_for_xlsx)]

# 将基因ID转换为家族名称
colnames(group_ratios_top)[-1] <- convert_gene_names(colnames(group_ratios_top)[-1])

# 确保类群按指定顺序排序
group_ratios_top <- group_ratios_top[order(group_ratios_top$Group), ]

# 转置数据框，使类群成为列，基因成为行
group_ratios_transposed <- as.data.frame(t(group_ratios_top[, -1]))
colnames(group_ratios_transposed) <- group_ratios_top$Group

# 添加基因家族名称作为行名
rownames(group_ratios_transposed) <- colnames(group_ratios_top)[-1]

# 将行名添加为单独的列
group_ratios_for_xlsx <- cbind(Gene_Family = rownames(group_ratios_transposed), 
                               group_ratios_transposed)

# 保存为xlsx文件
write_xlsx(list(Top_Genes_Group_Ratios = group_ratios_for_xlsx), 
           path = file.path(opt$output_dir, "top_genes_group_ratios.xlsx"))

cat("top", opt$top_n, "基因的类群平均成簇率已保存到:", 
    file.path(opt$output_dir, "top_genes_group_ratios.xlsx"), "\n")

cat("\n分析完成! 所有结果已保存到目录:", opt$output_dir, "\n")