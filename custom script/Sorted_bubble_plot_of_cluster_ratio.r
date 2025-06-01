#!/usr/bin/env Rscript
#########################################################################################################################
# 功能描述: 根据4类关键登陆性状的成簇率绘制气泡图，并根据连续成簇率为0的数量进行排序，便于比较不同家族的扩张情况
# 日期: 2025-05-02
# 作者：yyh
#########################################################################################################################

# 加载必要的R包
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

# 按性状类别定义颜色（新增list_6的颜色）
trait_colors <- c(
  "list_1" = "#fbb4ae",
  "list_2" = "#b3cde3",
  "list_3" = "#ccebc5",
  "list_6" = "#fed9a6"  # 为list_6添加新的颜色
)

# 类群展示顺序（从上到下）
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

# 提取性状和基因ID信息
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

# 5. 按性状过滤并排序，添加list_6性状
trait1_cols <- grep("^list_1_", gene_cols, value = TRUE)
trait2_cols <- grep("^list_2_", gene_cols, value = TRUE)
trait3_cols <- grep("^list_3_", gene_cols, value = TRUE)
trait6_cols <- grep("^list_6_", gene_cols, value = TRUE)  # 新增list_6

ratio_trait1 <- ratio_all[, c("Species", trait1_cols)]
ratio_trait2 <- ratio_all[, c("Species", trait2_cols)]
ratio_trait3 <- ratio_all[, c("Species", trait3_cols)]
ratio_trait6 <- ratio_all[, c("Species", trait6_cols)]  # 新增list_6
ratio_trait3_6 <- ratio_all[, c("Species", trait3_cols, trait6_cols)]  # 合并list_3和list_6

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

# 按类群聚合 - 分别处理不同性状类别，包括新增的list_6
group_ratio_trait1 <- calculate_group_ratio(ratio_trait1, group2spe)
group_ratio_trait2 <- calculate_group_ratio(ratio_trait2, group2spe)
group_ratio_trait3 <- calculate_group_ratio(ratio_trait3, group2spe)
group_ratio_trait6 <- calculate_group_ratio(ratio_trait6, group2spe)  # 新增list_6
group_ratio_trait3_6 <- calculate_group_ratio(ratio_trait3_6, group2spe)  # 合并list_3和list_6

# 按指定顺序排序类群
group_ratio_trait1$Group <- factor(group_ratio_trait1$Group, levels = group_order)
group_ratio_trait2$Group <- factor(group_ratio_trait2$Group, levels = group_order)
group_ratio_trait3$Group <- factor(group_ratio_trait3$Group, levels = group_order)
group_ratio_trait6$Group <- factor(group_ratio_trait6$Group, levels = group_order)  # 新增list_6
group_ratio_trait3_6$Group <- factor(group_ratio_trait3_6$Group, levels = group_order)  # 合并list_3和list_6

# 7. 将基因ID转换为家族名称
cat("将基因ID转换为家族名称...\n")
# 对类群聚合数据应用转换 - 分别处理不同性状类别，包括新增的list_6
group_ratio_trait1_family <- convert_col_names(group_ratio_trait1, gene_id_to_family, gene_order_map)
group_ratio_trait2_family <- convert_col_names(group_ratio_trait2, gene_id_to_family, gene_order_map)
group_ratio_trait3_family <- convert_col_names(group_ratio_trait3, gene_id_to_family, gene_order_map)
group_ratio_trait6_family <- convert_col_names(group_ratio_trait6, gene_id_to_family, gene_order_map)  # 新增list_6
group_ratio_trait3_6_family <- convert_col_names(group_ratio_trait3_6, gene_id_to_family, gene_order_map)  # 合并list_3和list_6

# 8. 保存数据为CSV - 更新为4类性状的数据
cat("保存数据为CSV...\n")
# 保存原始基因ID版本
write.csv(ratio_trait1, file.path(opt$output_dir, "gene_clustering_ratio_trait1_geneID.csv"), row.names = FALSE)
write.csv(ratio_trait2, file.path(opt$output_dir, "gene_clustering_ratio_trait2_geneID.csv"), row.names = FALSE)
write.csv(ratio_trait3, file.path(opt$output_dir, "gene_clustering_ratio_trait3_geneID.csv"), row.names = FALSE)
write.csv(ratio_trait6, file.path(opt$output_dir, "gene_clustering_ratio_trait6_geneID.csv"), row.names = FALSE)  # 新增list_6

# 保存类群聚合后的基因ID版本
write.csv(group_ratio_trait1, file.path(opt$output_dir, "group_clustering_ratio_trait1_geneID.csv"), row.names = FALSE)
write.csv(group_ratio_trait2, file.path(opt$output_dir, "group_clustering_ratio_trait2_geneID.csv"), row.names = FALSE)
write.csv(group_ratio_trait3, file.path(opt$output_dir, "group_clustering_ratio_trait3_geneID.csv"), row.names = FALSE)
write.csv(group_ratio_trait6, file.path(opt$output_dir, "group_clustering_ratio_trait6_geneID.csv"), row.names = FALSE)  # 新增list_6

# 保存转换为家族名称的版本
write.csv(group_ratio_trait1_family, file.path(opt$output_dir, "group_clustering_ratio_trait1_family.csv"), row.names = FALSE)
write.csv(group_ratio_trait2_family, file.path(opt$output_dir, "group_clustering_ratio_trait2_family.csv"), row.names = FALSE)
write.csv(group_ratio_trait3_family, file.path(opt$output_dir, "group_clustering_ratio_trait3_family.csv"), row.names = FALSE)
write.csv(group_ratio_trait6_family, file.path(opt$output_dir, "group_clustering_ratio_trait6_family.csv"), row.names = FALSE)  # 新增list_6

# 9. 保存为Excel文件（包含4类性状）
cat("保存数据为XLSX...\n")
write_xlsx(
  list(
    trait1_geneID      = ratio_trait1,
    trait2_geneID      = ratio_trait2,
    trait3_geneID      = ratio_trait3,
    trait6_geneID      = ratio_trait6,  # 新增list_6
    group_trait1_geneID = group_ratio_trait1,
    group_trait2_geneID = group_ratio_trait2,
    group_trait3_geneID = group_ratio_trait3,
    group_trait6_geneID = group_ratio_trait6  # 新增list_6
  ),
  path = file.path(opt$output_dir, "cluster_ratio_geneID.xlsx")
)

write_xlsx(
  list(
    trait1_family      = convert_col_names(ratio_trait1, gene_id_to_family, gene_order_map),
    trait2_family      = convert_col_names(ratio_trait2, gene_id_to_family, gene_order_map),
    trait3_family      = convert_col_names(ratio_trait3, gene_id_to_family, gene_order_map),
    trait6_family      = convert_col_names(ratio_trait6, gene_id_to_family, gene_order_map),  # 新增list_6
    group_trait1_family = group_ratio_trait1_family,
    group_trait2_family = group_ratio_trait2_family,
    group_trait3_family = group_ratio_trait3_family,
    group_trait6_family = group_ratio_trait6_family  # 新增list_6
  ),
  path = file.path(opt$output_dir, "cluster_ratio_family.xlsx")
)

# 10. 准备气泡图数据 - 计算基因家族的连续零值数量
cat("准备气泡图数据...\n")

# 计算连续零值数量函数
# 从Chlorophyta开始，计算每个基因家族有多少个连续的成簇率为0的类群
calculate_consecutive_zeros <- function(data, gene_family_col = "GeneFamily") {
  # 确保数据按照指定的类群顺序排序（从下到上：Chlorophyta到Angiosperms）
  data <- data %>% arrange(match(Group, rev(group_order)))
  
  # 获取所有基因家族
  gene_families <- unique(data[[gene_family_col]])
  
  # 存储每个基因家族的连续零值数量
  consecutive_zeros <- numeric(length(gene_families))
  names(consecutive_zeros) <- gene_families
  
  # 对每个基因家族计算连续零值的数量
  for (i in seq_along(gene_families)) {
    family <- gene_families[i]
    # 获取该基因家族的所有数据行
    family_data <- data[data[[gene_family_col]] == family, ]
    # 从Chlorophyta（底部）开始
    family_data <- family_data[order(match(family_data$Group, rev(group_order))), ]
    
    # 计算连续零值数量
    count <- 0
    for (j in 1:nrow(family_data)) {
      if (family_data$Ratio[j] == 0) {
        count <- count + 1
      } else {
        break  # 一旦遇到非零值，就停止计数
      }
    }
    
    consecutive_zeros[i] <- count
  }
  
  return(consecutive_zeros)
}

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
  
  # 创建排序信息
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

# 准备各个性状的气泡图数据，包括新增的list_6
bubble_data_trait1 <- prepare_bubble_data(group_ratio_trait1_family)
bubble_data_trait2 <- prepare_bubble_data(group_ratio_trait2_family)
bubble_data_trait3 <- prepare_bubble_data(group_ratio_trait3_family)
bubble_data_trait6 <- prepare_bubble_data(group_ratio_trait6_family)  # 新增list_6
bubble_data_trait3_6 <- prepare_bubble_data(group_ratio_trait3_6_family)  # 合并list_3和list_6

# 11. 按照连续零值数量分别排序各个性状的基因家族
cat("按照连续零值数量排序基因家族...\n")

# 分别计算各个性状基因家族的连续零值，包括新增的list_6
consecutive_zeros_trait1 <- calculate_consecutive_zeros(bubble_data_trait1)
consecutive_zeros_trait2 <- calculate_consecutive_zeros(bubble_data_trait2)
consecutive_zeros_trait3 <- calculate_consecutive_zeros(bubble_data_trait3)
consecutive_zeros_trait6 <- calculate_consecutive_zeros(bubble_data_trait6)  # 新增list_6

cat("Trait 1连续零值数量示例：\n")
print(head(consecutive_zeros_trait1))
cat("Trait 2连续零值数量示例：\n")
print(head(consecutive_zeros_trait2))
cat("Trait 3连续零值数量示例：\n")
print(head(consecutive_zeros_trait3))
cat("Trait 6连续零值数量示例：\n")  # 新增list_6
print(head(consecutive_zeros_trait6))

# 添加连续零值数量到数据集
bubble_data_trait1$ConsecutiveZeros <- sapply(bubble_data_trait1$GeneFamily, 
                                              function(x) consecutive_zeros_trait1[x])
bubble_data_trait2$ConsecutiveZeros <- sapply(bubble_data_trait2$GeneFamily, 
                                              function(x) consecutive_zeros_trait2[x])
bubble_data_trait3$ConsecutiveZeros <- sapply(bubble_data_trait3$GeneFamily, 
                                              function(x) consecutive_zeros_trait3[x])
bubble_data_trait6$ConsecutiveZeros <- sapply(bubble_data_trait6$GeneFamily, 
                                              function(x) consecutive_zeros_trait6[x])  # 新增list_6

# 获取排序后的基因家族
sorted_families_trait1 <- names(consecutive_zeros_trait1)[order(consecutive_zeros_trait1, decreasing = TRUE)]
sorted_families_trait2 <- names(consecutive_zeros_trait2)[order(consecutive_zeros_trait2, decreasing = TRUE)]
sorted_families_trait3 <- names(consecutive_zeros_trait3)[order(consecutive_zeros_trait3, decreasing = TRUE)]
sorted_families_trait6 <- names(consecutive_zeros_trait6)[order(consecutive_zeros_trait6, decreasing = TRUE)]  # 新增list_6

# 根据连续零值数量排序数据
bubble_data_trait1$GeneFamily <- factor(bubble_data_trait1$GeneFamily, 
                                        levels = sorted_families_trait1)
bubble_data_trait2$GeneFamily <- factor(bubble_data_trait2$GeneFamily, 
                                        levels = sorted_families_trait2)
bubble_data_trait3$GeneFamily <- factor(bubble_data_trait3$GeneFamily, 
                                        levels = sorted_families_trait3)
bubble_data_trait6$GeneFamily <- factor(bubble_data_trait6$GeneFamily, 
                                        levels = sorted_families_trait6)  # 新增list_6

# 12. 合并数据准备绘图
# 合并list_1和list_2的数据，保持各自的排序（用于子图a）
bubble_data_1_2 <- rbind(
  # 先放list_1的数据（排在左侧）
  bubble_data_trait1,
  # 再放list_2的数据（排在右侧）
  bubble_data_trait2
)

# 确保GeneFamily因子级别按照合并后的顺序
bubble_data_1_2$GeneFamily <- factor(bubble_data_1_2$GeneFamily, 
                                     levels = c(sorted_families_trait1, sorted_families_trait2))

# 合并list_3和list_6的数据，保持各自的排序（用于子图b）
bubble_data_3_6 <- rbind(
  # 先放list_3的数据（排在左侧）
  bubble_data_trait3,
  # 再放list_6的数据（排在右侧）
  bubble_data_trait6
)

# 确保GeneFamily因子级别按照合并后的顺序
bubble_data_3_6$GeneFamily <- factor(bubble_data_3_6$GeneFamily, 
                                     levels = c(sorted_families_trait3, sorted_families_trait6))

# 13. 创建气泡图 - 使用新的排序
cat("绘制气泡图...\n")
create_bubble_plot <- function(data, title = "", tag = "") {
  # 确保GeneFamily按新的排序方式排列，此时data$GeneFamily已经是一个factor，排序已由上一步完成
  
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
      legend.title = element_text(size = 7, angle = 0,face = "bold"),
      legend.title.position = "top",
      legend.direction = "vertical",
      legend.box            = "vertical",
      legend.box.just       = "center",                  
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
                         tag = "a"
)

p2 <- create_bubble_plot(bubble_data_3_6, 
                         title = "", 
                         tag = "b"
)

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

# 保存排序数据
cat("保存排序数据...\n")
sort_data_trait1 <- data.frame(
  GeneFamily = names(consecutive_zeros_trait1),
  ConsecutiveZeros = unname(consecutive_zeros_trait1),
  Trait = "list_1"
)
sort_data_trait1 <- sort_data_trait1[order(sort_data_trait1$ConsecutiveZeros, decreasing = TRUE), ]

sort_data_trait2 <- data.frame(
  GeneFamily = names(consecutive_zeros_trait2),
  ConsecutiveZeros = unname(consecutive_zeros_trait2),
  Trait = "list_2"
)
sort_data_trait2 <- sort_data_trait2[order(sort_data_trait2$ConsecutiveZeros, decreasing = TRUE), ]

sort_data_trait3 <- data.frame(
  GeneFamily = names(consecutive_zeros_trait3),
  ConsecutiveZeros = unname(consecutive_zeros_trait3),
  Trait = "list_3"
)
sort_data_trait3 <- sort_data_trait3[order(sort_data_trait3$ConsecutiveZeros, decreasing = TRUE), ]

sort_data_trait6 <- data.frame(
  GeneFamily = names(consecutive_zeros_trait6),
  ConsecutiveZeros = unname(consecutive_zeros_trait6),
  Trait = "list_6"
)
sort_data_trait6 <- sort_data_trait6[order(sort_data_trait6$ConsecutiveZeros, decreasing = TRUE), ]

write.csv(sort_data_trait1, file.path(opt$output_dir, "gene_family_sort_trait1.csv"), row.names = FALSE)
write.csv(sort_data_trait2, file.path(opt$output_dir, "gene_family_sort_trait2.csv"), row.names = FALSE)
write.csv(sort_data_trait3, file.path(opt$output_dir, "gene_family_sort_trait3.csv"), row.names = FALSE)
write.csv(sort_data_trait6, file.path(opt$output_dir, "gene_family_sort_trait6.csv"), row.names = FALSE)

# 另外保存一个合并后的排序信息文件
sort_data_all <- rbind(sort_data_trait1, sort_data_trait2, sort_data_trait3, sort_data_trait6)
write.csv(sort_data_all, file.path(opt$output_dir, "gene_family_sort_all.csv"), row.names = FALSE)

cat("处理完成! 所有结果已保存到", opt$output_dir, "目录\n")