#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(stringr)
library(optparse)
library(writexl)
# 定义命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="成簇基因分布的CSV文件路径", metavar="character"),
  make_option(c("-g", "--group"), type="character", default=NULL, 
              help="类群物种关系文件路径", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="基因家族名称映射文件路径", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="output", 
              help="输出目录", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 验证输入参数
if (is.null(opt$input) || is.null(opt$group) || is.null(opt$name)) {
  stop("提供输入文件路径")
}

# 创建输出目录
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# 按性状类别定义颜色
trait_colors <- c(
  "list_1" = "#fbb4ae",
  "list_2" = "#b3cde3",
  "list_3" = "#ccebc5"
)

# 类群展示顺序
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

# 1. 读取数据
cat("读取输入数据...\n")
gene_distribution <- read.csv(opt$input, check.names = FALSE)

# 读取类群物种映射
group2spe <- read.csv(opt$group, header = FALSE, 
                      col.names = c("Group", "Species"),
                      stringsAsFactors = FALSE)

# 读取基因家族映射
gene_family_map <- read.csv(opt$name, header = FALSE, 
                            col.names = c("Family", "GeneID"),
                            stringsAsFactors = FALSE)

# 创建基因ID到家族名称的映射字典
gene_id_to_family <- setNames(gene_family_map$Family, gene_family_map$GeneID)

# 打印映射信息用于调试
cat("读取到", nrow(gene_family_map), "个基因家族映射\n")
if(nrow(gene_family_map) > 0) {
  cat("示例映射：", gene_family_map$GeneID[1], "->", gene_family_map$Family[1], "\n")
}

# 2. 识别基因和成簇基因列
cat("识别基因列和计算成簇率...\n")
gene_cols <- grep("^list_[1-7]_.*[^_clustered]$", colnames(gene_distribution), value = TRUE)
clustered_cols <- grep("^list_[1-7]_.*_clustered$", colnames(gene_distribution), value = TRUE)

cat("找到", length(gene_cols), "个基因列和", length(clustered_cols), "个成簇基因列\n")

# 提取性状和基因ID信息（仅用于性状分类，不用于映射）
extract_trait <- function(gene_name) {
  pattern <- "^(list_[1-7])_.*$"
  matches <- str_match(gene_name, pattern)
  if (is.na(matches[1,1])) {
    return(NA)
  }
  return(matches[, 2])
}

# 创建一个全局顺序映射，记录每个基因ID的原始顺序
gene_order_map <- setNames(seq_along(gene_cols), gene_cols)

# 3. 计算成簇率
calculate_ratio <- function(df, gene_cols, clustered_cols) {
  ratio_data <- df[, "Species", drop = FALSE]
  
  for (i in seq_along(gene_cols)) {
    gene_col <- gene_cols[i]
    clustered_col <- clustered_cols[i]
    
    # 计算比例，处理分母为0的情况
    ratio_data[[gene_col]] <- ifelse(df[[gene_col]] > 0, 
                                     df[[clustered_col]] / df[[gene_col]], 
                                     0)
  }
  
  return(ratio_data)
}

# 计算全部七类性状的成簇率
ratio_all <- calculate_ratio(gene_distribution, gene_cols, clustered_cols)

# 4. 创建基因ID转换为家族名称的函数，同时保留原始顺序
convert_col_names <- function(ratio_df, id_to_name_map, order_map) {
  result_df <- ratio_df
  new_colnames <- colnames(ratio_df)
  
  # 创建一个映射记录原始列名和新列名的对应关系
  original_to_new <- list()
  
  # 跳过第一列（Species或Group）
  for (i in 2:length(new_colnames)) {
    original_colname <- new_colnames[i]
    
    # 直接使用完整列名查找匹配
    if (original_colname %in% names(id_to_name_map)) {
      # 使用家族名称替换列名
      new_colname <- id_to_name_map[original_colname]
      original_to_new[[original_colname]] <- new_colname
      new_colnames[i] <- new_colname
      cat("替换列名：", original_colname, " -> ", new_colname, "\n")
    } else {
      cat("警告：基因ID", original_colname, "没有对应的家族名称\n")
      original_to_new[[original_colname]] <- original_colname
    }
  }
  
  # 应用新的列名
  colnames(result_df) <- new_colnames
  
  # 保存原始顺序信息，为列名替换后的DF添加列属性
  attr(result_df, "original_to_new") <- original_to_new
  
  return(result_df)
}

# 5. 按性状过滤并排序
# 获取前三类性状的列
trait1_3_cols <- grep("^list_[1-3]_", gene_cols, value = TRUE)
ratio_trait1_3 <- ratio_all[, c("Species", trait1_3_cols)]

# 6. 按类群聚合数据
cat("按类群聚合数据...\n")
calculate_group_ratio <- function(ratio_data, group2spe) {
  # 创建group到species的映射
  species_to_group <- setNames(group2spe$Group, group2spe$Species)
  
  # 添加Group列
  ratio_data$Group <- species_to_group[ratio_data$Species]
  
  # 移除NA的Group（可能有物种不在映射中）
  ratio_data <- ratio_data[!is.na(ratio_data$Group), ]
  
  # 按Group聚合并计算平均值
  group_ratio <- ratio_data %>%
    group_by(Group) %>%
    summarise(across(-Species, ~mean(.x, na.rm = TRUE)))
  
  return(group_ratio)
}

# 按类群聚合
group_ratio_all <- calculate_group_ratio(ratio_all, group2spe)
group_ratio_1_3 <- calculate_group_ratio(ratio_trait1_3, group2spe)

# 按指定顺序排序类群
group_ratio_all$Group <- factor(group_ratio_all$Group, levels = group_order)
group_ratio_1_3$Group <- factor(group_ratio_1_3$Group, levels = group_order)

# 7. 将基因ID转换为家族名称
cat("将基因ID转换为家族名称...\n")
# 对类群聚合数据应用转换
group_ratio_all_family <- convert_col_names(group_ratio_all, gene_id_to_family, gene_order_map)
group_ratio_1_3_family <- convert_col_names(group_ratio_1_3, gene_id_to_family, gene_order_map)

# 对物种级数据应用转换
ratio_all_family <- convert_col_names(ratio_all, gene_id_to_family, gene_order_map)
ratio_trait1_3_family <- convert_col_names(ratio_trait1_3, gene_id_to_family, gene_order_map)

# 8. 保存数据为CSV - (基因ID版本 & 家族名称版本)
cat("保存数据为CSV...\n")
write.csv(ratio_all, file.path(opt$output_dir, "gene_clustering_ratio_all_geneID.csv"), row.names = FALSE)
write.csv(ratio_trait1_3, file.path(opt$output_dir, "gene_clustering_ratio_1_3_geneID.csv"), row.names = FALSE)
write.csv(group_ratio_all, file.path(opt$output_dir, "group_clustering_ratio_all_geneID.csv"), row.names = FALSE)
write.csv(group_ratio_1_3, file.path(opt$output_dir, "group_clustering_ratio_1_3_geneID.csv"), row.names = FALSE)

write.csv(ratio_all_family, file.path(opt$output_dir, "gene_clustering_ratio_all_family.csv"), row.names = FALSE)
write.csv(ratio_trait1_3_family, file.path(opt$output_dir, "gene_clustering_ratio_1_3_family.csv"), row.names = FALSE)
write.csv(group_ratio_all_family, file.path(opt$output_dir, "group_clustering_ratio_all_family.csv"), row.names = FALSE)
write.csv(group_ratio_1_3_family, file.path(opt$output_dir, "group_clustering_ratio_1_3_family.csv"), row.names = FALSE)

# 9. 保存为两个 Excel 文件（基因ID版本 & 家族名称版本）
cat("保存数据为 XLSX...\n")
write_xlsx(
  list(
    all_geneID        = ratio_all,
    trait1_3_geneID   = ratio_trait1_3,
    group_all_geneID  = group_ratio_all,
    group1_3_geneID   = group_ratio_1_3
  ),
  path = file.path(opt$output_dir, "cluster_ratio_geneID.xlsx")
)

write_xlsx(
  list(
    all_family        = ratio_all_family,
    trait1_3_family   = ratio_trait1_3_family,
    group_all_family  = group_ratio_all_family,
    group1_3_family   = group_ratio_1_3_family
  ),
  path = file.path(opt$output_dir, "cluster_ratio_family.xlsx")
)

# 10. 准备气泡图数据 - 保持原始基因顺序
cat("准备气泡图数据...\n")
prepare_bubble_data <- function(group_ratio_df) {
  # 获取原始列顺序信息
  original_to_new <- attr(group_ratio_df, "original_to_new")
  if (is.null(original_to_new)) {
    # 如果没有顺序映射，创建一个基于当前列名的映射
    original_to_new <- setNames(colnames(group_ratio_df)[-1], colnames(group_ratio_df)[-1])
  }
  
  # 创建反向映射（新列名到原始列名）
  new_to_original <- list()
  for (orig in names(original_to_new)) {
    new_name <- original_to_new[[orig]]
    new_to_original[[new_name]] <- orig
  }
  
  # 创建排序信息 - 基于原始基因ID顺序
  gene_order <- list()
  for (new_name in colnames(group_ratio_df)[-1]) {
    # 找到对应的原始列名
    original_name <- if (new_name %in% names(new_to_original)) 
      new_to_original[[new_name]] else new_name
    
    # 获取原始顺序编号
    if (original_name %in% names(gene_order_map)) {
      gene_order[[new_name]] <- gene_order_map[[original_name]]
    } else {
      # 如果找不到顺序，给一个较大的值（放在最后）
      gene_order[[new_name]] <- 1000 + length(gene_order)
    }
  }
  
  # 从宽格式转为长格式
  long_data <- group_ratio_df %>%
    pivot_longer(cols = -Group, names_to = "GeneFamily", values_to = "Ratio")
  
  # 添加顺序信息供后续排序
  long_data$OrderIndex <- sapply(long_data$GeneFamily, function(x) {
    if (x %in% names(gene_order)) {
      return(gene_order[[x]])
    }
    return(1000) # 找不到顺序信息的默认放在最后
  })
  
  # 添加性状信息
  long_data$Trait <- sapply(long_data$GeneFamily, function(x) {
    # 尝试从对应的原始基因ID中提取性状
    if (x %in% names(new_to_original)) {
      original_name <- new_to_original[[x]]
      trait <- extract_trait(original_name)
      if (!is.na(trait)) return(trait)
    }
    
    # 直接尝试从名称中提取性状
    trait <- extract_trait(x)
    if (!is.na(trait)) return(trait)
    
    # 反向查找 - 查找是否是基因家族名称
    for (gene_id in names(gene_id_to_family)) {
      if (gene_id_to_family[gene_id] == x) {
        trait <- extract_trait(gene_id)
        if (!is.na(trait)) return(trait)
      }
    }
    
    return("unknown")
  })
  
  return(long_data)
}

# 准备所有数据用于气泡图
bubble_data_all <- prepare_bubble_data(group_ratio_all_family)

# 调试信息
cat("气泡图数据中的唯一基因家族名称：\n")
print(head(unique(bubble_data_all$GeneFamily), 10))
cat("气泡图数据中的唯一性状类别：\n")
print(unique(bubble_data_all$Trait))

# 准备前三类性状的子集并按原始顺序排序
bubble_data_1_2 <- bubble_data_all %>% 
  filter(Trait %in% c("list_1", "list_2")) %>%
  arrange(OrderIndex)

bubble_data_3 <- bubble_data_all %>% 
  filter(Trait == "list_3") %>%
  arrange(OrderIndex)

# 11. 创建气泡图 - 保留原始基因顺序
cat("绘制气泡图...\n")
create_bubble_plot <- function(data, title = "", tag = "") {
  # 确保GeneFamily按原始顺序排列
  data$GeneFamily <- factor(data$GeneFamily, levels = unique(data$GeneFamily))
  
  # 绘制气泡图
  p <- ggplot(data, aes(x = GeneFamily, y = Group)) +
    geom_point(aes(size = Ratio * 100, fill = Trait), shape = 21, 
               color = "black", alpha = 0.8, stroke = 0.5) +
    scale_size_area(max_size = 7, name = "Clustering Rate (%)") +
    scale_fill_manual(values = trait_colors,
                      name = "Trait Category",
                      guide="none") +
    labs(title = title, x = "", y = "", tag = tag) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 270, vjust = 0, hjust = 1, size = 6, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7, angle = 90),
      legend.title.position = "left",
      legend.direction = "vertical",
      plot.tag = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 12)
    ) +
    scale_x_discrete(position = "top") +  # x轴标签放在顶部
    scale_y_discrete(limits = rev(group_order))  # 确保类群按指定顺序显示
  
  return(p)
}

# 创建子图
p1 <- create_bubble_plot(bubble_data_1_2,
                         title = "",
                         tag = "a")

p2 <- create_bubble_plot(bubble_data_3,
                         title = "",
                         tag = "b")

# 合并子图
cat("生成最终图表...\n")
combined_plot <- plot_grid(
  p1, p2, 
  nrow = 2,
  rel_heights = c(1.1, 1),  # 调整高度比例
  align = 'v'
)

# 保存图表
ggsave(file.path(opt$output_dir, "cluster_ratio_bubble_plot.pdf"), 
       combined_plot, width = 8.3, height = 6.5, units = "in")
ggsave(file.path(opt$output_dir, "cluster_ratio_bubble_plot.png"), 
       combined_plot, width = 8.3, height = 6.5, units = "in", dpi = 300)

cat("处理完成! 所有结果已保存到", opt$output_dir, "目录\n")