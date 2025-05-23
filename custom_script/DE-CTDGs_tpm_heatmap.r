# 使用ComplexHeatmap绘制两个物种的TPM表达量聚类热图
# Drawing TPM expression clustered heatmaps for two species using ComplexHeatmap

# 加载必要的库 / Load necessary libraries
library(ComplexHeatmap)
library(grid)
library(circlize)  # 用于颜色映射 / For color mapping

# 设置输出文件名和格式 / Set output filename and format
output_file <- "TPM_expression_heatmaps_complex.pdf"
# A4尺寸 (in inches) / A4 size (in inches)
pdf_width <- 8.3
pdf_height <- 8

# 检查文件是否存在 / Check if files exist
check_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("错误：文件不存在 -", file_path, "\nError: File does not exist -", file_path))
  }
}

check_file("Aty_TPM_matrix.csv")
check_file("Ppa_TPM_matrix.csv")

# 读取CSV文件 / Read CSV files
aty_data <- read.csv("Aty_TPM_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
ppa_data <- read.csv("Ppa_TPM_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 显示数据结构信息 / Display data structure information
cat("Aty TPM矩阵结构 / Aty TPM matrix structure:\n")
cat("行数 (Rows):", nrow(aty_data), "- 基因/基因簇数量 / Number of genes/clusters\n")
cat("列数 (Columns):", ncol(aty_data), "- 处理条件数量 / Number of treatment conditions\n")
cat("处理条件名称 / Treatment condition names:\n")
print(colnames(aty_data))

cat("\nPpa TPM矩阵结构 / Ppa TPM matrix structure:\n")
cat("行数 (Rows):", nrow(ppa_data), "- 基因/基因簇数量 / Number of genes/clusters\n")
cat("列数 (Columns):", ncol(ppa_data), "- 处理条件数量 / Number of treatment conditions\n")
cat("处理条件名称 / Treatment condition names:\n")
print(colnames(ppa_data))

# 数据预处理：log2变换并添加伪计数 / Data preprocessing: log2 transform with pseudocount
# 直接应用log2(TPM+1)转换，不进行Z-score标准化 / Directly apply log2(TPM+1) transformation without Z-score normalization
aty_log2 <- log2(aty_data + 1)
ppa_log2 <- log2(ppa_data + 1)

# 设置颜色映射 - 从白色到蓝色的渐变 / Set color mapping - gradient from white to blue
col_fun <- colorRamp2(c(0, max(c(max(aty_log2), max(ppa_log2)))/2, max(c(max(aty_log2), max(ppa_log2)))), 
                      c("white", "#6baed6", "#2b83ba"))

# 设置PDF输出 / Set PDF output
pdf(output_file, width = pdf_width, height = pdf_height)

# 创建两个热图的列表 / Create a list of two heatmaps
ht_list <- list()

# 第一个热图 (Aty) / First heatmap (Aty)
ht1 <- Heatmap(as.matrix(aty_log2),
               name = "log2(TPM+1)",      # 图例名称 / Legend name
               col = col_fun,            # 颜色设置 / Color settings
               show_row_names = FALSE,    # 不显示行名 / Don't show row names
               border = FALSE,            # 无边框 / No border
               column_names_gp = gpar(fontsize = 12),  # 列名字体大小 / Column name font size
               cluster_rows = TRUE,       # 对行进行聚类 / Cluster rows
               cluster_columns = FALSE,   # 不对列进行聚类 / Don't cluster columns
               column_names_side = "top", # 列名在顶部 / Column names on top
               column_names_rot = 0,     # 列名旋转角度 / Column name rotation
               show_heatmap_legend = TRUE, # 显示图例 / Show legend
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10, fontface = "bold"),
                 labels_gp = gpar(fontsize = 8)
               )
)

# 第二个热图 (Ppa) / Second heatmap (Ppa)
ht2 <- Heatmap(as.matrix(ppa_log2),
               name = "log2(TPM+1)",      # 图例名称 / Legend name
               col = col_fun,            # 颜色设置 / Color settings
               show_row_names = FALSE,    # 不显示行名 / Don't show row names
               border = FALSE,            # 无边框 / No border
               column_names_gp = gpar(fontsize = 12),  # 列名字体大小 / Column name font size
               cluster_rows = TRUE,       # 对行进行聚类 / Cluster rows
               cluster_columns = FALSE,   # 不对列进行聚类 / Don't cluster columns
               column_names_side = "top", # 列名在顶部 / Column names on top
               column_names_rot = 0,     # 列名旋转角度 / Column name rotation
               show_heatmap_legend = TRUE, # 显示图例 / Show legend
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10, fontface = "bold"),
                 labels_gp = gpar(fontsize = 8)
               )
)

# 设置页面布局并绘制热图 / Set page layout and draw heatmaps
# 创建一个1列2行的布局 / Create a 1-column 2-row layout
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 1, heights = unit(c(1, 1), "null"))))

# 绘制第一个热图 (Aty) / Draw first heatmap (Aty)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht1, newpage = FALSE)
# 添加左上角标签 'a' / Add label 'a' in the top-left corner
grid.text("a", x = 0.02, y = 0.98, just = c("left", "top"), 
          gp = gpar(fontsize = 16, fontface = "bold"))
popViewport()

# 绘制第二个热图 (Ppa) / Draw second heatmap (Ppa)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht2, newpage = FALSE)
# 添加左上角标签 'b' / Add label 'b' in the top-left corner
grid.text("b", x = 0.02, y = 0.98, just = c("left", "top"), 
          gp = gpar(fontsize = 16, fontface = "bold"))
popViewport()

# 关闭PDF设备 / Close PDF device
dev.off()

cat(paste("\n\nTPM表达量热图已保存至", output_file, "\nTPM expression heatmaps saved to", output_file, "\n"))