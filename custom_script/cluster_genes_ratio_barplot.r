# 植物基因组成簇基因分析与可视化脚本
# 从CSV文件中读取220个物种的基因组数据，并创建垂直排列的柱状图


# 加载必要的库
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})

# 设置输出选项
output_file <- "clustered_genes_ratio.pdf"
create_segments <- FALSE                      # 是否创建分段图
segment_size <- 44                           # 每个分段的物种数量（分成5段，每段约44个物种）
color_scheme <- "magma"                    # 可选: "viridis", "plasma", "inferno", "magma", "cividis"
label_decimal_places <- 2                    # 标签中小数位数

# 函数：读取并处理数据
process_data <- function(file_path) {
  # 读取CSV文件
  cat("读取数据文件:", file_path, "\n")
  data <- tryCatch({
    read_csv(file_path, col_types = cols(.default = "c"))  # 首先全部作为字符串读取
  }, error = function(e) {
    stop("无法读取CSV文件: ", e$message)
  })
  
  # 检查列名
  expected_cols <- c("species", "number of cluestered genes", "number of genes", "percentage")
  actual_cols <- names(data)
  
  # 如果列名不匹配，尝试修复
  if (!all(expected_cols %in% actual_cols)) {
    cat("警告: CSV文件列名与预期不符\n")
    cat("预期列名: ", paste(expected_cols, collapse = ", "), "\n")
    cat("实际列名: ", paste(actual_cols, collapse = ", "), "\n")
    
    if (ncol(data) >= 4) {
      cat("尝试按位置分配列名\n")
      names(data)[1:4] <- expected_cols
    } else {
      stop("CSV文件列数不足")
    }
  }
  
  # 创建清理后的数据框
  data_clean <- data %>%
    rename(
      species = `species`,
      clustered_genes = `number of cluestered genes`,
      total_genes = `number of genes`,
      percentage = `percentage`
    )
  
  # 确保数据类型正确
  data_clean <- data_clean %>%
    mutate(
      clustered_genes = as.numeric(clustered_genes),
      total_genes = as.numeric(total_genes),
      percentage = as.numeric(percentage)
    )
  
  # 检查并处理缺失值
  na_rows <- which(is.na(data_clean$percentage))
  if (length(na_rows) > 0) {
    cat("发现", length(na_rows), "行的百分比值缺失，尝试重新计算\n")
    
    # 尝试从原始数据重新计算
    data_clean <- data_clean %>%
      mutate(percentage = ifelse(is.na(percentage) & !is.na(clustered_genes) & !is.na(total_genes),
                                 (clustered_genes / total_genes) * 100, 
                                 percentage))
    
    # 检查是否仍有NA值
    remaining_na <- which(is.na(data_clean$percentage))
    if (length(remaining_na) > 0) {
      warning("以下行的百分比值仍无法计算: ", 
              paste(remaining_na, collapse = ", "))
    }
  }
  
  # 保持物种的原始顺序
  data_clean$species <- factor(data_clean$species, levels = data_clean$species)
  
  # 返回处理后的数据
  cat("数据处理完成: 包含", nrow(data_clean), "个物种\n")
  return(data_clean)
}

# 函数：创建颜色方案
create_color_palette <- function(n, scheme = "viridis") {
  if (scheme == "viridis") {
    return(viridis(n, option = "D"))
  } else if (scheme == "plasma") {
    return(viridis(n, option = "A"))
  } else if (scheme == "inferno") {
    return(viridis(n, option = "B"))
  } else if (scheme == "magma") {
    return(viridis(n, option = "C"))
  } else if (scheme == "cividis") {
    return(viridis(n, option = "E"))
  } else {
    # 默认回退到viridis
    warning("未知的颜色方案: ", scheme, "，使用默认的viridis")
    return(viridis(n, option = "D"))
  }
}

# 函数：创建主图
create_main_plot <- function(data, color_palette, decimal_places = 2) {
  # 创建格式化字符串
  format_string <- paste0("%.", decimal_places, "f%%")
  
  # 创建垂直柱状图
  p <- ggplot(data, aes(x = percentage, y = species, fill = species)) +
    geom_col() +  # 创建柱状图
    scale_fill_manual(values = color_palette) +  # 应用自定义颜色
    geom_text(aes(label = sprintf(format_string, percentage)), 
              hjust = -0.5, size = 1.5) +  # 在柱状图右侧添加百分比标签
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(), 
      panel.grid.minor.y = element_blank(),
      axis.line.x = element_line(linewidth = .3,color = "black"),
      axis.ticks.x = element_line(linewidth = .3,color = "black"),
      axis.ticks.y = element_line(linewidth = .3,color = "black"),
      axis.text.y = element_text(size = 3.5,color = "black",angle = 0,hjust = 1), 
      plot.margin = margin(0, 0, 0, 40), 
      axis.title.y.left = element_blank(),
      axis.title.x.bottom = element_blank(),
      legend.position = "none"  
    ) +
    labs(
      y = "",
      title = "",
      subtitle = ""
    ) +
    # 确保柱状图紧贴坐标轴
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(data$percentage, na.rm = TRUE) * 1.2))
  
  return(p)
}

# 函数：创建分段图
create_segment_plots <- function(data, color_palette, segment_size, decimal_places = 2) {
  # 计算分段数量
  n_species <- nrow(data)
  num_segments <- ceiling(n_species / segment_size)
  
  # 创建格式化字符串
  format_string <- paste0("%.", decimal_places, "f%%")
  
  # 创建每个分段的图
  for (i in 1:num_segments) {
    # 计算此分段的起始和结束索引
    start_idx <- (i - 1) * segment_size + 1
    end_idx <- min(i * segment_size, n_species)
    
    # 提取此分段的数据
    segment_data <- data[start_idx:end_idx, ]
    
    # 创建此分段的图
    p_segment <- ggplot(segment_data, aes(x = percentage, y = species, fill = species)) +
      geom_col() +
      scale_fill_manual(values = color_palette[start_idx:end_idx]) +
      geom_text(aes(label = sprintf(format_string, percentage)), 
                hjust = -0.1, size = 3) +  # 可以使用更大的字体
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size = 8),  # 可以使用更大的字体
        plot.margin = margin(10, 30, 10, 10),
        legend.position = "none"
      ) +
      labs(
        x = "成簇基因百分比 (%)",
        title = paste0("各物种基因成簇率 - 部分 ", i, "/", num_segments),
        subtitle = paste0("物种 ", start_idx, " 到 ", end_idx, " (共", n_species, "个)")
      ) +
      scale_x_continuous(expand = c(0, 0), limits = c(0, max(data$percentage, na.rm = TRUE) * 1.2))
    
    # 保存此分段
    output_file <- paste0("clustered_genes_ratio_part_", i, ".pdf")
    cat("保存分段图 ", i, "/", num_segments, " 到文件: ", output_file, "\n")
    
    ggsave(output_file, p_segment, width = 8.3, height = 11.7, units = "in", dpi = 300)
  }
  
  cat("已创建", num_segments, "个分段图\n")
}

# 主程序
main <- function() {
  # 设置开始时间
  start_time <- Sys.time()
  cat("开始处理: ", format(start_time), "\n")
  
  # 读取并处理数据
  data <- process_data("cluster_genes_ratio.csv")
  
  # 创建颜色方案
  cat("使用颜色方案: ", color_scheme, "\n")
  color_palette <- create_color_palette(nrow(data), scheme = color_scheme)
  
  # 创建主图
  cat("创建主图...\n")
  main_plot <- create_main_plot(data, color_palette, decimal_places = label_decimal_places)
  
  # 保存主图
  cat("保存主图到文件: ", output_file, "\n")
  ggsave(output_file, main_plot, width = 8.3, height = 11.7, units = "in", 
         dpi = 300, limitsize = FALSE)
  
  # 创建分段图（如果需要）
  if (create_segments) {
    cat("创建分段图...\n")
    create_segment_plots(data, color_palette, segment_size, decimal_places = label_decimal_places)
  }
  
  # 计算运行时间
  end_time <- Sys.time()
  run_time <- difftime(end_time, start_time, units = "secs")
  cat("处理完成，用时: ", round(run_time, 2), " 秒\n")
}

# 执行主程序
main()