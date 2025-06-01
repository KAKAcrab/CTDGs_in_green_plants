#!/usr/bin/env Rscript

# DE-CTDGs成员表达量可视化与分析脚本
# 作者: YYH
# 日期: 2025-04-14
#
# 功能：
# 1. 读取并处理每种胁迫条件下的*_member_expression.csv文件
# 2. 为每种胁迫条件创建配对箱线图
# 3. 使用配对Wilcoxon检验对major和minor之间的表达量差异进行显著性检验
# 4. 将不同胁迫条件的结果用子图展示
# 5. 在各个子图左上角添加小写字母标注
#
# 使用方法:
# Rscript plot_member_expression.R --input_dir /path/to/input_directory --output_dir /path/to/output_directory

# 加载必要的R包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(argparse)
  library(gridExtra)
  library(ggpubr)
  library(grid)
})

# 解析命令行参数
parser <- ArgumentParser(description="DE-CTDGs成员表达量可视化脚本")
parser$add_argument("--input_dir", dest="input_dir", required=TRUE, help="包含*_member_expression.csv文件的输入目录")
parser$add_argument("--output_dir", dest="output_dir", required=TRUE, help="输出图像文件的目录")
parser$add_argument("--file_pattern", dest="file_pattern", default="*_member_expression.csv", help="要处理的文件模式, 默认为'*_member_expression.csv'")
parser$add_argument("--width", dest="width", type="double", default=8.3, help="输出图像宽度（英寸）")
parser$add_argument("--height", dest="height", type="double", default=4, help="输出图像高度（英寸）")
parser$add_argument("--dpi", dest="dpi", type="integer", default=300, help="输出图像分辨率，默认为300")
parser$add_argument("--format", dest="format", default="pdf", help="输出图像格式, 默认为'pdf'")
parser$add_argument("--y_label", dest="y_label", default="Expression (TPM)", help="Y轴标签，默认为'Expression (TPM)'")
parser$add_argument("--adjust_method", dest="adjust_method", default="BH", help="多重检验校正方法, 默认为'BH'(Benjamini-Hochberg)")

args <- parser$parse_args()

# 创建输出目录（如果不存在）
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}

# 获取输入文件列表
file_pattern <- file.path(args$input_dir, args$file_pattern)
input_files <- list.files(path = args$input_dir, pattern = basename(args$file_pattern), full.names = TRUE)

if (length(input_files) == 0) {
  stop("未找到匹配的输入文件: ", file_pattern)
}

# 读取并处理所有输入文件
cat("读取", length(input_files), "个输入文件...\n")

# 函数：读取并转换单个文件
process_file <- function(file_path) {
  # 从文件名中提取条件名（去掉_member_expression.csv部分）
  condition <- gsub("_member_expression\\.csv$", "", basename(file_path))
  
  # 读取CSV文件
  data <- read_csv(file_path, show_col_types = FALSE)
  
  # 检查数据格式
  required_cols <- c("cluster_id", "major_gene", "minor_gene", "major_tpm", "minor_tpm")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    warning(paste("文件缺少必要的列:", paste(missing_cols, collapse=", ")))
    return(NULL)
  }
  
  # 转换为长格式，便于绘图
  data_long <- rbind(
    data %>% select(cluster_id, gene = major_gene, tpm = major_tpm) %>% mutate(type = "Major"),
    data %>% select(cluster_id, gene = minor_gene, tpm = minor_tpm) %>% mutate(type = "Minor")
  )
  
  # 添加条件信息
  data_long$condition <- condition
  
  return(data_long)
}

# 处理所有文件并合并数据
all_data <- map_dfr(input_files, process_file)

if (is.null(all_data) || nrow(all_data) == 0) {
  stop("未能成功处理任何输入文件")
}

# 确保type列是因子类型，并设置level顺序
all_data$type <- factor(all_data$type, levels = c("Major", "Minor"))

# 获取所有条件
conditions <- unique(all_data$condition)
cat("找到", length(conditions), "个条件:", paste(conditions, collapse=", "), "\n")

# 使用配对Wilcoxon检验并获取显著性标记
get_significance <- function(condition_data) {
  # 转换为宽格式数据以配对观测
  wide_data <- condition_data %>%
    pivot_wider(id_cols = cluster_id, names_from = type, values_from = tpm)
  
  # 检查配对数据是否完整
  if (any(is.na(wide_data$Major)) || any(is.na(wide_data$Minor))) {
    message("数据中存在NA值，将被移除用于配对检验")
    # 移除包含NA值的行
    wide_data <- wide_data %>% filter(!is.na(Major), !is.na(Minor))
  }
  
  # 确保有足够的配对样本
  if (nrow(wide_data) < 2) {
    warning("配对样本数量不足，无法执行Wilcoxon检验")
    return(list(p_value = NA, sig_mark = "NA", test_name = "无法检验"))
  }
  
  # 运行配对Wilcoxon符号秩检验
  tryCatch({
    wilcox_result <- wilcox.test(wide_data$Major, wide_data$Minor, paired = TRUE)
    test_name <- "Wilcoxon配对检验"
    p_value <- wilcox_result$p.value
  }, error = function(e) {
    # 如果Wilcoxon检验失败（例如，所有差异为0），回退到配对t检验
    message("Wilcoxon检验失败，回退到配对t检验: ", e$message)
    t_test_result <- t.test(wide_data$Major, wide_data$Minor, paired = TRUE)
    return(list(p_value = t_test_result$p.value, 
                sig_mark = get_sig_mark(t_test_result$p.value),
                test_name = "配对t检验(回退)"))
  })
  
  # 获取显著性标记
  sig_mark <- get_sig_mark(p_value)
  
  list(p_value = p_value, sig_mark = sig_mark, test_name = test_name)
}

# 辅助函数：根据p值获取显著性标记
get_sig_mark <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# 计算每个条件的显著性
sig_results <- all_data %>%
  group_by(condition) %>%
  group_map(~get_significance(.x), .keep = TRUE) %>%
  setNames(conditions)

# 如果有多个条件，进行多重检验校正
if (length(conditions) > 1) {
  # 收集所有p值
  p_values <- sapply(sig_results, function(x) x$p_value)
  
  # 使用指定方法进行校正
  adjusted_p <- p.adjust(p_values, method = args$adjust_method)
  
  # 更新sig_results中的p值和显著性标记
  for (i in seq_along(sig_results)) {
    sig_results[[i]]$original_p_value <- sig_results[[i]]$p_value
    sig_results[[i]]$p_value <- adjusted_p[i]
    sig_results[[i]]$sig_mark <- get_sig_mark(adjusted_p[i])
    sig_results[[i]]$adjusted <- TRUE
  }
  
  cat("已使用", args$adjust_method, "方法对", length(p_values), "个p值进行多重检验校正\n")
}

# 为每个条件创建一个子图，并添加小写字母标签
plot_condition <- function(condition, index) {
  # 提取该条件的数据
  condition_data <- all_data %>% filter(condition == !!condition)
  
  # 获取显著性结果
  sig_result <- sig_results[[condition]]
  
  # 计算最大值，用于放置显著性标记
  y_max <- max(condition_data$tpm) * 1.1
  
  # 美化条件名称（替换下划线为空格，首字母大写）
  plot_title <- gsub("_", " ", condition)
  plot_title <- tools::toTitleCase(plot_title)
  
  # 创建小写字母标签（从a开始）
  panel_label <- letters[index]
  
  # 绘图
  # 计算每个type的箱线图统计信息（四分位数和 whisker范围）
  pts_filtered <- condition_data %>%
    group_by(type) %>%
    mutate(
      Q1 = quantile(tpm, 0.25, na.rm = TRUE),
      Q3 = quantile(tpm, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower = Q1 - 1.5 * IQR,
      upper = Q3 + 1.5 * IQR
    ) %>%
    ungroup() %>%
    # 只保留落在 whisker 范围内的点
    filter(tpm >= lower, tpm <= upper)
  
  # 确定y轴范围使用过滤后的数据
  y_min <- min(pts_filtered$tpm, na.rm = TRUE)
  y_max <- max(pts_filtered$tpm, na.rm = TRUE)
  
  # 箱线图\用完整数据,不显示 outlier，散点与连线使用过滤后的数据
  p <- ggplot(condition_data, aes(x = type, y = tpm)) +
    # 绘制中空箱线图（边框颜色由 type 决定，内部填充设为NA）
    geom_boxplot(aes(color = type), fill = NA,
                 width = 0.5, alpha = 0.7, linewidth = 0.4,
                 outlier.shape = NA) +
    # 仅显示箱线图内的点,无需抖动,保持点在同一x位置
    geom_point(data = pts_filtered, aes(color = type),
               position = position_identity(), size = 0.2, alpha = 0.7) +
    # 连接过滤后点的连线，确保只连接留下来的点
    geom_line(data = pts_filtered, aes(group = cluster_id),
              color = "grey", alpha = 0.5, linewidth = 0.2) +
    # 添加显著性标记
    annotate("text", x = 1.5, y = y_max, label = sig_result$sig_mark, size = 3) +
    # 添加小写字母标签
    annotate("text", x = 0.6, y = y_max * 0.95, label = panel_label,
             size = 4, fontface = "bold", hjust = 1) +
    # 坐标轴和标题
    labs(
      title = plot_title,
      x = "",
      y = args$y_label
    ) +
    # 设置颜色
    scale_color_manual(values = c("Major" = "#e41a1c", "Minor" = "#56B4E9")) +
    # 使用主题，移除不必要的背景与网格
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 8, color = "black"),
      axis.title.y = element_text(size = 10, color = "black"),
      axis.line.y = element_line(linewidth = 0.3, color = "black"),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 3, r = 3, b = 3, l = 3, unit = "pt")
    ) +
    # 设置 y 轴范围，正好覆盖过滤后点的范围
    scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0))
  
  return(p)
}

# 创建所有条件的子图
plot_list <- map2(conditions, seq_along(conditions), plot_condition)

# 计算行列数（尽量使布局接近正方形）
n_plots <- length(conditions)
n_cols <- ceiling(sqrt(n_plots))
n_rows <- ceiling(n_plots / n_cols)

# 设置更紧凑的布局
plot_layout <- grid.layout(
  nrow = n_rows,
  ncol = n_cols,
  widths = unit(rep(1, n_cols), "null"),
  heights = unit(rep(1, n_rows), "null"),
  respect = FALSE
)

# 合并所有子图，使用更紧凑的布局
if (length(plot_list) > 1) {
  # 创建一个新的绘图设备
  combined_plot <- arrangeGrob(
    grobs = plot_list, 
    ncol = n_cols,
    top = textGrob(
      "Expression Differences Between Major and Minor members in DE-CTDGs of P.patens",
      gp = gpar(fontsize = 8, fontface = "bold")
    ),
    padding = unit(0.2, "line")  # 减小子图之间的间距
  )
} else {
  combined_plot <- plot_list[[1]]
}

# 保存图像
output_file <- file.path(args$output_dir, paste0("member_expression_comparison.", args$format))
ggsave(
  filename = output_file,
  plot = combined_plot,
  width = args$width,
  height = args$height,
  dpi = args$dpi
)

cat("图像已保存到:", output_file, "\n")

# 保存统计检验结果
test_results <- data.frame(
  condition = names(sig_results),
  test_name = sapply(sig_results, function(x) x$test_name),
  p_value = sapply(sig_results, function(x) x$p_value),
  significance = sapply(sig_results, function(x) x$sig_mark)
)

# 添加原始p值（如果进行了多重检验校正）
if (length(conditions) > 1) {
  test_results$original_p_value <- sapply(sig_results, function(x) 
    if("original_p_value" %in% names(x)) x$original_p_value else x$p_value)
  test_results$adjustment_method <- args$adjust_method
}

test_output <- file.path(args$output_dir, "statistical_test_results.csv")
write_csv(test_results, test_output)

cat("统计检验结果已保存到:", test_output, "\n")

# 创建有关ratio分布的箱线图
ratio_data <- input_files %>%
  map_dfr(function(file_path) {
    condition <- gsub("_member_expression\\.csv$", "", basename(file_path))
    data <- read_csv(file_path, show_col_types = FALSE)
    data$condition <- condition
    return(data)
  })

# 计算比率并绘制箱线图
if ("ratio" %in% colnames(ratio_data)) {
  p_ratio <- ggplot(ratio_data, aes(x = condition, y = ratio)) +
    geom_boxplot(fill = "#56B4E9", alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(
      title = "Distribution of Expression Ratio (Major/Minor) Across Conditions",
      x = "",
      y = "Expression Ratio (Major/Minor)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8, color = "black"),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # 保存比率分布图
  ratio_output_file <- file.path(args$output_dir, paste0("ratio_distribution.", args$format))
  ggsave(
    filename = ratio_output_file,
    plot = p_ratio,
    width = args$width * 0.8,
    height = args$height * 0.6,
    dpi = args$dpi
  )
  
  cat("表达比率分布图已保存到:", ratio_output_file, "\n")
}