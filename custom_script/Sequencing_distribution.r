#!/usr/bin/env Rscript

# 载入所需包
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(optparse)

# 定义命令行参数
option_list <- list(
  make_option(c("-t", "--taxon"), type="character", default=NULL, 
              help="物种类群映射文件路径 (abbr_taxon)", metavar="character"),
  make_option(c("-a", "--assembly"), type="character", default=NULL, 
              help="测序类型文件路径 (assembly_type)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="sequencing_distribution.pdf", 
              help="输出PDF文件名 [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 验证输入参数
if (is.null(opt$taxon) || is.null(opt$assembly)) {
  stop("请提供物种类群映射文件和测序类型文件路径")
}

# 1. 读取数据
cat("读取输入数据...\n")
taxon_data <- read.csv(opt$taxon, header = FALSE, 
                       col.names = c("Species", "Group"),
                       stringsAsFactors = FALSE)
assembly_data <- read.csv(opt$assembly, header = FALSE, 
                          col.names = c("Species", "Assembly_Type"),
                          stringsAsFactors = FALSE)

# 2. 数据合并和预处理
cat("处理数据...\n")
combined_data <- merge(taxon_data, assembly_data, by = "Species")

# 定义类群顺序（按进化关系）
group_order <- c("Chlorophyta", "Streptophyte_algae", "Zygnematophyceae", 
                 "Bryophytes", "Lycophytes", "Ferns", 
                 "Gymnosperms", "Angiosperms")

# 定义测序类型顺序和颜色
assembly_order <- c("TGS", "NGS+TGS", "FGS+NGS+TGS", "NGS", "FGS+NGS", "FGS")

# 颜色映射：按读长技术能力分组
assembly_colors <- c(
  "TGS" = "#3182bd",           # 深蓝：纯长读长
  "NGS+TGS" = "#6baed6",       # 中蓝：长读长混合
  "FGS+NGS+TGS" = "#bdd7e7",   # 浅蓝：全技术混合
  "NGS" = "#eff3ff",           # 浅蓝：纯短读长
  "FGS+NGS" = "#fee5d9",       # 灰色：一代+短读长
  "FGS" = "#fcae91"            # 深灰：纯一代测序
)

# 确保因子顺序
combined_data$Group <- factor(combined_data$Group, levels = group_order)
combined_data$Assembly_Type <- factor(combined_data$Assembly_Type, levels = assembly_order)

# 3. 计算总体分布（饼图数据）
cat("计算总体分布...\n")
overall_dist <- combined_data %>%
  count(Assembly_Type) %>%
  mutate(prop = n / sum(n) * 100,
         label = paste0(Assembly_Type, "\n(", round(prop, 1), "%)"))

# 4. 计算各类群分布（柱状图数据）
cat("计算各类群分布...\n")
group_dist <- combined_data %>%
  group_by(Group, Assembly_Type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(prop = count / sum(count) * 100)

# 5. 创建饼图
create_pie_chart <- function(data) {
  # 计算标签位置
  data <- data %>%
    arrange(Assembly_Type) %>%
    mutate(
      cumsum = cumsum(prop),
      pos = cumsum - prop/2
    )
  
  ggplot(data, aes(x = "", y = prop, fill = Assembly_Type)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = assembly_colors) +
    geom_text(aes(y = pos, label = paste0(round(prop, 1), "%")), 
              color = "white", fontface = "bold", size = 3) +
    labs(title = "Overall Distribution") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "none"
    )
}

# 6. 创建堆叠柱状图
create_group_chart <- function(data) {
  ggplot(data, aes(x = Group, y = prop, fill = Assembly_Type)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = assembly_colors,
      name  = "Assembly Type",
      guide = guide_legend(
        direction      = "vertical",  # 图例项垂直排列
        title.position = "top",       # 标题放在图例项上方
        title.hjust    = 0.5          # 标题水平居中
      )
    ) +
    scale_x_discrete(limits = group_order) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      title = "Distribution by Taxonomic Groups",
      x     = "",
      y     = "Proportion (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.x     = element_text(angle = 30, hjust = 1, size = 8, color = "black"),
      axis.text.y     = element_text(size = 9, color = "black"),
      axis.title.y    = element_text(size = 10),
      axis.line.y     = element_line(linewidth = .3, color = "black"),
      axis.line.x     = element_blank(),
      axis.ticks      = element_line(linewidth = .3, color = "black"),
      legend.position = "left",      # 图例放右侧
      legend.direction= "vertical",   # 竖直排列
      legend.title    = element_text(size = 10),
      legend.text     = element_text(size = 9),
      panel.grid      = element_blank(),
      panel.border    = element_blank()
    )
}

# 7. 生成图表
cat("生成图表...\n")
pie_plot <- create_pie_chart(overall_dist)
group_plot <- create_group_chart(group_dist)

# 8. 组合图表
combined_plot <- plot_grid(
  pie_plot, group_plot,
  nrow = 1,
  rel_widths = c(1, 1.5),
  align = 'h'
)

# 9. 创建综合统计表格
cat("生成统计表格...\n")

# 9.1 创建详细的数量和比例表格
create_summary_table <- function(combined_data, group_order, assembly_order) {
  # 计算各类群各测序类型的详细统计
  detailed_stats <- combined_data %>%
    group_by(Group, Assembly_Type) %>%
    summarise(count = n(), .groups = "drop") %>%
    # 确保所有类群-测序类型组合都存在
    tidyr::complete(Group = factor(group_order, levels = group_order),
                    Assembly_Type = factor(assembly_order, levels = assembly_order),
                    fill = list(count = 0)) %>%
    # 计算各类群比例
    group_by(Group) %>%
    mutate(group_total = sum(count),
           group_prop = ifelse(group_total > 0, 
                               round(count / group_total * 100, 1), 
                               0)) %>%
    ungroup()
  
  # 计算总体统计
  overall_stats <- combined_data %>%
    count(Assembly_Type) %>%
    mutate(overall_total = sum(n),
           overall_prop = round(n / overall_total * 100, 1)) %>%
    # 确保所有测序类型都存在
    tidyr::complete(Assembly_Type = factor(assembly_order, levels = assembly_order),
                    fill = list(n = 0)) %>%
    mutate(overall_total = sum(n, na.rm = TRUE),
           overall_prop = round(n / overall_total * 100, 1))
  
  # 创建宽格式表格（各类群为行，测序类型为列）
  count_table <- detailed_stats %>%
    select(Group, Assembly_Type, count) %>%
    tidyr::pivot_wider(names_from = Assembly_Type, values_from = count, values_fill = 0) %>%
    mutate(Group_Total = rowSums(select(., -Group)))
  
  prop_table <- detailed_stats %>%
    select(Group, Assembly_Type, group_prop) %>%
    tidyr::pivot_wider(names_from = Assembly_Type, values_from = group_prop, values_fill = 0) %>%
    mutate(Group_Total = 100.0)  # 每个类群的比例总和应为100%
  
  # 添加总体行
  overall_count <- overall_stats %>%
    select(Assembly_Type, n) %>%
    tidyr::pivot_wider(names_from = Assembly_Type, values_from = n, values_fill = 0) %>%
    mutate(Group_Total = rowSums(across(everything()))) %>%
    mutate(Group = "Overall", .before = 1)
  
  overall_prop <- overall_stats %>%
    select(Assembly_Type, overall_prop) %>%
    tidyr::pivot_wider(names_from = Assembly_Type, values_from = overall_prop, values_fill = 0) %>%
    mutate(Group_Total = 100.0) %>%
    mutate(Group = "Overall", .before = 1)
  
  # 合并数量表格
  count_table_final <- bind_rows(count_table, overall_count)
  
  # 合并比例表格
  prop_table_final <- bind_rows(prop_table, overall_prop)
  
  return(list(
    count_table = count_table_final,
    prop_table = prop_table_final,
    detailed_stats = detailed_stats,
    overall_stats = overall_stats
  ))
}

# 生成统计表格
summary_tables <- create_summary_table(combined_data, group_order, assembly_order)

# 9.2 保存统计表格为CSV文件
output_prefix <- tools::file_path_sans_ext(opt$output)

# 保存数量表格
write.csv(summary_tables$count_table, 
          paste0(output_prefix, "_count_table.csv"), 
          row.names = FALSE)

# 保存比例表格
write.csv(summary_tables$prop_table, 
          paste0(output_prefix, "_proportion_table.csv"), 
          row.names = FALSE)

# 保存详细统计表格
detailed_table <- summary_tables$detailed_stats %>%
  select(Group, Assembly_Type, count, group_prop) %>%
  arrange(factor(Group, levels = group_order), 
          factor(Assembly_Type, levels = assembly_order))

write.csv(detailed_table, 
          paste0(output_prefix, "_detailed_stats.csv"), 
          row.names = FALSE)

# 10. 保存为A4尺寸PDF
cat("保存PDF文件...\n")
ggsave(opt$output, combined_plot, 
       width = 8.3, height = 5, units = "in",
       device = "pdf")

cat("完成! 输出文件:", opt$output, "\n")

# 11. 输出数据统计摘要
cat("\n=== 数据统计摘要 ===\n")
cat("总物种数:", nrow(combined_data), "\n")
cat("包含的类群:", length(unique(combined_data$Group)), "\n")
cat("测序技术类型:", length(unique(combined_data$Assembly_Type)), "\n")

cat("\n总体分布:\n")
print(overall_dist %>% select(Assembly_Type, n, prop))

cat("\n各类群样本数:\n")
group_summary <- combined_data %>% 
  count(Group) %>% 
  arrange(factor(Group, levels = group_order))
print(group_summary)

cat("\n=== 统计表格已生成 ===\n")
cat("数量表格:", paste0(output_prefix, "_count_table.csv"), "\n")
cat("比例表格:", paste0(output_prefix, "_proportion_table.csv"), "\n") 
cat("详细统计:", paste0(output_prefix, "_detailed_stats.csv"), "\n")

# 12. 显示数量表格预览
cat("\n数量表格预览:\n")
print(summary_tables$count_table)

cat("\n比例表格预览:\n")
print(summary_tables$prop_table)