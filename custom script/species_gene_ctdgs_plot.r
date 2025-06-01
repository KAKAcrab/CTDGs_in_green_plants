#!/usr/bin/env Rscript

# 载入所需包 / Load required packages
library(dplyr)
library(ggplot2)
library(optparse)

# 定义命令行参数 / Define command line arguments
option_list <- list(
  make_option(c("-c", "--ctdgs"), type="character", default=NULL, 
              help="CTDGs数量文件路径 (CTDGs_number_for_220_spe)", metavar="character"),
  make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="基因数量文件路径 (gene_number_for_220_spe)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="species_gene_ctdgs_plot.pdf", 
              help="输出PDF文件名 [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 验证输入参数 / Validate input parameters
if (is.null(opt$ctdgs) || is.null(opt$genes)) {
  stop("请提供CTDGs数量文件和基因数量文件路径 / Please provide CTDGs and gene count file paths")
}

# 1. 读取数据 / Read data
cat("读取输入数据... / Reading input data...\n")
ctdgs_data <- read.csv(opt$ctdgs, header = FALSE, 
                       col.names = c("Species", "CTDGs_Count"),
                       stringsAsFactors = FALSE)
genes_data <- read.csv(opt$genes, header = FALSE, 
                       col.names = c("Species", "Gene_Count"),
                       stringsAsFactors = FALSE)

# 2. 数据合并和预处理 / Merge and preprocess data
cat("处理数据... / Processing data...\n")
combined_data <- merge(ctdgs_data, genes_data, by = "Species")

# 检查数据完整性 / Check data integrity
if (nrow(combined_data) == 0) {
  stop("合并后数据为空，请检查两个文件中的物种名称是否匹配 / No data after merging, please check species names")
}

# 3. 按基因数量降序排序 / Sort by gene count in descending order
combined_data <- combined_data %>%
  arrange(desc(Gene_Count)) %>%
  mutate(Species = factor(Species, levels = Species))  # 固定因子顺序 / Fix factor order

# 4. 数据统计信息 / Data statistics
cat("\n=== 数据统计信息 / Data Statistics ===\n")
cat("物种数量 / Number of species:", nrow(combined_data), "\n")
cat("基因数量范围 / Gene count range:", min(combined_data$Gene_Count), "-", max(combined_data$Gene_Count), "\n")
cat("CTDGs数量范围 / CTDGs count range:", min(combined_data$CTDGs_Count), "-", max(combined_data$CTDGs_Count), "\n")

# 5. 创建组合图表 / Create combined plot
cat("生成图表... / Generating plots...\n")

# 计算次坐标轴比例因子 / Calculate secondary axis scaling factor
# 使得CTDGs的最大值约等于基因数量最大值的80%，保持视觉平衡
scale_factor <- max(combined_data$Gene_Count) * 0.8 / max(combined_data$CTDGs_Count)

create_combined_plot <- function(data) {
  ggplot(data, aes(x = Species)) +
    # 绘制基因数量柱状图 / Draw gene count bar chart
    geom_bar(aes(y = Gene_Count), 
             stat = "identity", 
             fill = "#3498db", 
             alpha = 0.7, 
             width = 1) +
    
    # 绘制CTDGs折线图 / Draw CTDGs line chart
    geom_line(aes(y = CTDGs_Count * scale_factor, group = 1), 
              color = "#e74c3c", 
              linewidth = .5) +
    
    # 绘制CTDGs散点 / Draw CTDGs points
    geom_point(aes(y = CTDGs_Count * scale_factor), 
               color = "#e74c3c", 
               size = .8) +
    
    # 设置主y轴（基因数量）/ Set primary y-axis (gene count)
    scale_y_continuous(
      name = "Gene Count",
      expand = c(0, 0),  # y轴从0开始 / y-axis starts from 0
      limits = c(0, max(data$Gene_Count) * 1.1),  # 留出10%上边距 / 10% top margin
      labels = function(x) format(x, scientific = FALSE, big.mark = ","),
      
      # 设置次y轴（CTDGs数量）/ Set secondary y-axis (CTDGs count)
      sec.axis = sec_axis(~ . / scale_factor, 
                         name = "CTDGs Count",
                         labels = function(x) format(round(x), scientific = FALSE, big.mark = ","))
    ) +
    
    # x轴标签 / x-axis labels
    scale_x_discrete(name = "Species") +
    
    # 应用简洁主题 / Apply minimal theme
    theme_minimal() +
    theme(
      # 文本设置 / Text settings
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      axis.title.y.right = element_text(size = 12, margin = margin(l = 10)),
      
      # 坐标轴设置 / Axis settings
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      
      # 网格设置 / Grid settings
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      # 背景设置 / Background settings
      panel.background = element_blank(),
      plot.background = element_blank(),
      
      # 边距设置 / Margin settings
      plot.margin = margin(20, 20, 20, 20)
    ) +
    
    # 添加图例说明 / Add legend annotation
    annotation_custom(
      grob = grid::textGrob(
        label = paste0("Blue bars: Gene count\n",
                      "Red line: CTDGs count"),
        x = 0.02, y = 0.98,
        hjust = 0, vjust = 1,
        gp = grid::gpar(fontsize = 10, col = "black")
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
}

# 生成图表 / Generate plot
final_plot <- create_combined_plot(combined_data)

# 6. 保存为横版A4尺寸PDF / Save as landscape A4 PDF
cat("保存PDF文件... / Saving PDF file...\n")
ggsave(opt$output, final_plot, 
       width = 11.7, height = 8.3, units = "in",  # 横版A4尺寸 / Landscape A4 size
       device = "pdf")

cat("完成! 输出文件: / Completed! Output file:", opt$output, "\n")

# 7. 输出详细统计摘要 / Output detailed statistics summary
cat("\n=== 详细统计摘要 / Detailed Statistics Summary ===\n")

# 基因数量统计 / Gene count statistics
cat("基因数量统计 / Gene Count Statistics:\n")
cat("  平均值 / Mean:", round(mean(combined_data$Gene_Count)), "\n")
cat("  中位数 / Median:", round(median(combined_data$Gene_Count)), "\n")
cat("  标准差 / Standard deviation:", round(sd(combined_data$Gene_Count)), "\n")

# CTDGs数量统计 / CTDGs count statistics
cat("CTDGs数量统计 / CTDGs Count Statistics:\n")
cat("  平均值 / Mean:", round(mean(combined_data$CTDGs_Count)), "\n")
cat("  中位数 / Median:", round(median(combined_data$CTDGs_Count)), "\n")
cat("  标准差 / Standard deviation:", round(sd(combined_data$CTDGs_Count)), "\n")

# 相关性分析 / Correlation analysis
correlation <- cor(combined_data$Gene_Count, combined_data$CTDGs_Count)
cat("基因数量与CTDGs数量的相关系数 / Correlation between gene and CTDGs counts:", round(correlation, 3), "\n")

# 前5和后5的物种 / Top 5 and bottom 5 species
cat("\n基因数量最多的5个物种 / Top 5 species by gene count:\n")
top5 <- head(combined_data, 5)
for(i in 1:nrow(top5)) {
  cat("  ", i, ".", top5$Species[i], "- 基因/Genes:", top5$Gene_Count[i], 
      "CTDGs:", top5$CTDGs_Count[i], "\n")
}

cat("\n基因数量最少的5个物种 / Bottom 5 species by gene count:\n")
bottom5 <- tail(combined_data, 5)
for(i in 1:nrow(bottom5)) {
  cat("  ", i, ".", bottom5$Species[i], "- 基因/Genes:", bottom5$Gene_Count[i], 
      "CTDGs:", bottom5$CTDGs_Count[i], "\n")
}