#!/usr/bin/env Rscript

#############################################################
# CTDG-FoldChange Visualization Script
# Author: Claude
# Date: 2025-05-20
#
# 创建CTDG-FoldChange结果的可视化。
# 该版本专注于跨条件的组合图。
#
# 主要功能:
# 1. 组合密度图: 比较DE-CTDGs和随机CTDGs之间的fold change分布
# 2. 组合箱形图: 表达分布的视觉比较
#
# 使用方法:
#   Rscript visualize_ctdg_foldchange.R \
#     --data_dir /path/to/analysis_output \
#     --output_dir /path/to/output \
#     [--size_groups "member=2,3<member<5,6<member<10,member>10"] \
#     [--color_palette "viridis"] \
#     [--width 8.3] \
#     [--height 4]
#############################################################

# 加载所需包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)
  library(viridis)  # 用于颜色方案
})

# 解析命令行参数
parse_arguments <- function() {
  option_list <- list(
    make_option("--data_dir", type="character", default=NULL, 
                help="包含Python分析脚本输出的数据目录"),
    make_option("--output_dir", type="character", default="./plots", 
                help="生成图表的输出目录"),
    make_option("--size_groups", type="character", default="member=2,3<member<5,6<member<10,member>10",
                help="要处理的大小组，逗号分隔"),
    make_option("--color_palette", type="character", default="viridis",
                help="颜色方案，选项：viridis, brewer, grey"),
    make_option("--width", type="integer", default=8.3,
                help="图表宽度（英寸）"),
    make_option("--height", type="integer", default=4,
                help="图表高度（英寸）")
  )
  
  parser <- OptionParser(option_list=option_list,
                         description="可视化CTDG-FoldChange分析结果")
  
  args <- parse_args(parser)
  
  # 检查必需参数
  if (is.null(args$data_dir)) {
    stop("必须提供数据目录 (--data_dir)")
  }
  
  # 从conditions.txt读取条件列表（如果存在）
  conditions_file <- file.path(args$data_dir, "conditions.txt")
  if (file.exists(conditions_file)) {
    args$conditions <- readLines(conditions_file)
    cat("从 conditions.txt 文件读取条件列表:", paste(args$conditions, collapse=", "), "\n")
  } else {
    # 从visualization_data目录推断条件
    viz_dir <- file.path(args$data_dir, "visualization_data")
    if (dir.exists(viz_dir)) {
      args$conditions <- list.dirs(viz_dir, full.names=FALSE, recursive=FALSE)
      if (length(args$conditions) > 0) {
        cat("从目录推断条件:", paste(args$conditions, collapse=", "), "\n")
      } else {
        stop("无法从目录推断条件，且未找到conditions.txt文件")
      }
    } else {
      stop("无法从目录推断条件，且未找到conditions.txt文件")
    }
  }
  
  # 创建输出目录
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive=TRUE)
  }
  
  # 解析大小组
  args$size_groups <- strsplit(args$size_groups, ",")[[1]]
  
  return(args)
}

# 加载可视化数据
load_visualization_data <- function(data_dir, condition, verbose=TRUE) {
  if (verbose) cat("加载条件", condition, "的数据\n")
  
  viz_dir <- file.path(data_dir, "visualization_data", condition)
  
  if (!dir.exists(viz_dir)) {
    warning("目录不存在: ", viz_dir)
    return(NULL)
  }
  
  # 加载密度图数据
  density_file <- file.path(viz_dir, "density_data.csv")
  if (file.exists(density_file)) {
    density_data <- read_csv(density_file, show_col_types=FALSE)
    if (verbose) cat("  已加载密度数据:", nrow(density_data), "行\n")
  } else {
    warning("密度数据文件未找到: ", density_file)
    density_data <- NULL
  }
  
  # 加载箱形图数据
  boxplot_file <- file.path(viz_dir, "boxplot_data.csv")
  if (file.exists(boxplot_file)) {
    boxplot_data <- read_csv(boxplot_file, show_col_types=FALSE)
    if (verbose) cat("  已加载箱形图数据:", nrow(boxplot_data), "行\n")
  } else {
    warning("箱形图数据文件未找到: ", boxplot_file)
    boxplot_data <- NULL
  }
  
  return(list(
    density_data = density_data,
    boxplot_data = boxplot_data
  ))
}

# 创建组合密度图
create_combined_density_plot <- function(all_density_data, size_groups, conditions, color_palette) {
  # 检查是否有数据
  if (length(all_density_data) == 0) {
    warning("没有可用的密度数据来创建组合图")
    return(NULL)
  }
  
  # 将所有数据合并为一个数据框
  combined_data <- data.frame()
  
  for (condition in names(all_density_data)) {
    data <- all_density_data[[condition]]
    if (!is.null(data) && nrow(data) > 0) {
      # 筛选指定大小组的数据
      data <- data %>% filter(size_group %in% size_groups)
      
      # 添加条件列
      data$condition <- condition
      
      # 添加到组合数据
      combined_data <- rbind(combined_data, data)
    }
  }
  
  if (nrow(combined_data) == 0) {
    warning("筛选和合并后没有有效数据")
    return(NULL)
  }
  
  # 设置大小组的正确顺序
  size_group_order <- c("member=2", "3<member<5", "6<member<10", "member>10")
  combined_data$size_group <- factor(combined_data$size_group, levels=size_group_order)
  
  # 设置条件的正确顺序
  combined_data$condition <- factor(combined_data$condition, levels=conditions)
  
  # 创建大小组的标签
  size_group_labels <- c(
    "member=2" = "Members = 2",
    "3<member<5" = "Members 3-5",
    "6<member<10" = "Members 6-10",
    "member>10" = "Members >10"
  )
  
  # 创建组合图
  p <- ggplot(combined_data, aes(x=log2_fold_change, fill=group)) +
    geom_density(alpha=0.7, size=.1) +
    facet_grid(condition ~ size_group, 
               scales="free_y", 
               labeller=labeller(size_group=as_labeller(size_group_labels))) +
    scale_fill_manual(values=c(
      "DE-CTDG" = "#e41a1c",    
      "Random"  = "#377eb8"), 
      name="Groups",
      labels=c("DE-CTDG" = "DE-CTDGs", "Random" = "Random CTDGs")) +
    theme_minimal() +
    labs(title="",
         x="Log2 Fold Change",
         y="Density") +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size=12),
      axis.text = element_text(size=8, color = "black"),
      axis.line.x = element_line(linewidth = .3, color = "black"),
      axis.line.y = element_line(linewidth = .3, color = "black"),
      legend.title = element_text(size=12),
      legend.text = element_text(size=10),
      strip.text = element_text(size=11, face="bold"),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    geom_vline(xintercept=0, linetype="dashed", color="black", linewidth=.3)
  
  return(p)
}

# 创建组合箱形图
create_combined_boxplot <- function(all_boxplot_data, size_groups, conditions, color_palette) {
  # 检查是否有数据
  if (length(all_boxplot_data) == 0) {
    warning("没有可用的箱形图数据来创建组合图")
    return(NULL)
  }
  
  # 将所有数据合并为一个数据框
  combined_data <- data.frame()
  
  for (condition in names(all_boxplot_data)) {
    data <- all_boxplot_data[[condition]]
    if (!is.null(data) && nrow(data) > 0) {
      # 筛选指定大小组的数据
      data <- data %>% filter(size_group %in% size_groups)
      
      # 添加条件列
      data$condition <- condition
      
      # 添加到组合数据
      combined_data <- rbind(combined_data, data)
    }
  }
  
  if (nrow(combined_data) == 0) {
    warning("筛选和合并后没有有效数据")
    return(NULL)
  }
  
  # 设置大小组的正确顺序
  size_group_order <- c("member=2", "3<member<5", "6<member<10", "member>10")
  combined_data$size_group <- factor(combined_data$size_group, levels=size_group_order)
  
  # 设置条件的正确顺序
  combined_data$condition <- factor(combined_data$condition, levels=conditions)
  
  # 创建大小组的标签
  size_group_labels <- c(
    "member=2" = "Members = 2",
    "3<member<5" = "Members 3-5",
    "6<member<10" = "Members 6-10",
    "member>10" = "Members >10"
  )
  
  # 创建组合箱形图
  p <- ggplot(combined_data, aes(x=group, y=log2_fold_change, color=group)) +
    geom_boxplot(alpha=0.8, outlier.size=.3, width=0.6, fill=NA, linewidth=.3) +
    facet_grid(condition ~ size_group, 
               scales="free_y", 
               labeller=labeller(size_group=as_labeller(size_group_labels))) +
    scale_color_manual(
      values = c(
        "DE-CTDG" = "#e41a1c",    
        "Random"  = "#377eb8"),
      name   = NULL,
      guide  = "none"
    ) +
    theme_minimal() +
    labs(title="",
         x="Groups",
         y="Log2 Fold Change") +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title = element_text(size=12),
      axis.text.x = element_text(size=10, color = "black"),
      axis.text.y = element_text(size=7, color = "black"),
      axis.line.x = element_line(linewidth = .3, color = "black"),
      axis.line.y = element_line(linewidth = .3, color = "black"),
      legend.title = element_text(size=12),
      legend.text = element_text(size=10),
      strip.text = element_text(size=11, face="bold"),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=.3)
  
  return(p)
}

# 保存图表 (仅PDF)
save_plot <- function(plot, filename, width, height) {
  # 仅保存PDF版本
  tryCatch({
    pdf(filename, width=width, height=height)
    print(plot)
    dev.off()
    
    cat("图表已保存至:", filename, "\n")
  }, error = function(e) {
    warning("保存图表出错: ", e$message)
  })
}

# 主函数
main <- function() {
  # 解析命令行参数
  args <- parse_arguments()
  
  cat("数据目录:", args$data_dir, "\n")
  cat("输出目录:", args$output_dir, "\n")
  cat("处理条件:", paste(args$conditions, collapse=", "), "\n")
  cat("处理大小组:", paste(args$size_groups, collapse=", "), "\n")
  
  # 创建输出目录
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 处理每个条件，仅收集数据
  all_density_data <- list()
  all_boxplot_data <- list()
  
  for (condition in args$conditions) {
    cat("\n处理条件:", condition, "\n")
    
    # 加载数据
    viz_data <- load_visualization_data(args$data_dir, condition)
    
    if (is.null(viz_data)) {
      warning("跳过条件 ", condition, ": 无法加载数据")
      next
    }
    
    # 存储组合图的数据
    if (!is.null(viz_data$density_data) && nrow(viz_data$density_data) > 0) {
      all_density_data[[condition]] <- viz_data$density_data
    }
    
    if (!is.null(viz_data$boxplot_data) && nrow(viz_data$boxplot_data) > 0) {
      all_boxplot_data[[condition]] <- viz_data$boxplot_data
    }
  }
  
  # 创建组合密度图
  if (length(all_density_data) > 0) {
    tryCatch({
      combined_density <- create_combined_density_plot(
        all_density_data,
        args$size_groups,
        args$conditions,
        args$color_palette
      )
      
      if (!is.null(combined_density)) {
        # 保存组合密度图
        combined_density_file <- file.path(args$output_dir, "combined_density_plot.pdf")
        save_plot(combined_density, combined_density_file, args$width, args$height)
      }
    }, error = function(e) {
      warning("创建组合密度图出错: ", e$message)
    })
  }
  
  # 创建组合箱形图
  if (length(all_boxplot_data) > 0) {
    tryCatch({
      combined_boxplot <- create_combined_boxplot(
        all_boxplot_data,
        args$size_groups,
        args$conditions,
        args$color_palette
      )
      
      if (!is.null(combined_boxplot)) {
        # 保存组合箱形图
        combined_boxplot_file <- file.path(args$output_dir, "combined_boxplot_plot.pdf")
        save_plot(combined_boxplot, combined_boxplot_file, args$width, args$height)
      }
    }, error = function(e) {
      warning("创建组合箱形图出错: ", e$message)
    })
  }
  
  cat("\n所有条件处理完成\n")
}

# 执行主函数
main()