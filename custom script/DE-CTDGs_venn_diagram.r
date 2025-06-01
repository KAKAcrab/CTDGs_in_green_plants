# VennDiagram绘制DE-CTDGs的脚本（改进版）
# Improved script for creating Venn diagrams of DE-CTDGs

# 加载必要的库 / Load necessary libraries
library(VennDiagram)
library(gridExtra)
library(grid)

# 抑制VennDiagram包的信息消息（更简洁的输出）/ Suppress info messages from VennDiagram (cleaner output)
futile.logger::flog.threshold(futile.logger::ERROR)

# 设置输出文件名和格式 / Set output filename and format
output_file <- "DE_CTDGs_venn_diagram_improved.pdf"
# A4尺寸 (in inches) / A4 size (in inches)
pdf_width <- 8.3
pdf_height <- 8

# 检查文件是否存在 / Check if files exist
check_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("错误：文件不存在 -", file_path, "\nError: File does not exist -", file_path))
  }
}

check_file("Aty_DE_CTDGs_matrix.csv")
check_file("Ppa_DE_CTDGs_matrix.csv")

# 读取CSV文件 / Read CSV files
aty_data <- read.csv("Aty_DE_CTDGs_matrix.csv", header = TRUE, stringsAsFactors = FALSE)
ppa_data <- read.csv("Ppa_DE_CTDGs_matrix.csv", header = TRUE, stringsAsFactors = FALSE)

# 获取列名 / Get column names
aty_cols <- colnames(aty_data)
ppa_cols <- colnames(ppa_data)

# 获取基因簇ID的列名 / Get column name for gene cluster ID
cluster_id_col_aty <- aty_cols[1]
cluster_id_col_ppa <- ppa_cols[1]

# 获取处理条件的列名 / Get column names for treatment conditions
treatment_cols_aty <- aty_cols[3:length(aty_cols)]  # 跳过第一列(ID)和第二列(成员数量) / Skip first column (ID) and second column (member count)
treatment_cols_ppa <- ppa_cols[3:length(ppa_cols)]

cat("Aty处理条件 / Aty treatment conditions:", length(treatment_cols_aty), "个条件 / conditions\n")
print(treatment_cols_aty)
cat("\nPpa处理条件 / Ppa treatment conditions:", length(treatment_cols_ppa), "个条件 / conditions\n")
print(treatment_cols_ppa)

# 提取每个处理条件的基因集 / Extract gene sets for each treatment condition
# 对于Aty文件 / For Aty file
aty_sets <- list()
for (col in treatment_cols_aty) {
  aty_sets[[col]] <- aty_data[[cluster_id_col_aty]][aty_data[[col]] == 1]
  cat(paste("\n处理条件 / Treatment condition:", col, "- 基因数量 / Number of genes:", length(aty_sets[[col]])))
}

# 对于Ppa文件 / For Ppa file
cat("\n")
ppa_sets <- list()
for (col in treatment_cols_ppa) {
  ppa_sets[[col]] <- ppa_data[[cluster_id_col_ppa]][ppa_data[[col]] == 1]
  cat(paste("\n处理条件 / Treatment condition:", col, "- 基因数量 / Number of genes:", length(ppa_sets[[col]])))
}

# 设置颜色 / Set colors - 使用ColorBrewer配色方案 / Using ColorBrewer color scheme
colors_aty <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6")  # Set1
colors_ppa <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4")

# 设置PDF输出 / Set PDF output
pdf(output_file, width = pdf_width, height = pdf_height)

# 创建两个子图 / Create two subplots
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))

# 设置Venn图参数 / Set Venn diagram parameters
venn_params <- list(
  alpha = 0.5,
  col = "black",
  cex = .8,
  fontface = "bold",
  cat.cex = 1.0,
  cat.fontface = "plain",
  cat.dist = 0.20,  # 增加标签距离以避免遮蔽 / Increase label distance to avoid overlap
  margin = 0.1,
  sub.cex = 1.2,
  lwd = 0.1  # 添加线条宽度参数，使描边更细 / Add line width parameter to make borders thinner
)

# 第一个子图 (Aty) / First subplot (Aty)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
# 针对5个集合的优化布局 / Optimized layout for 5 sets
venn_aty <- venn.diagram(
  x = aty_sets,
  filename = NULL,
  fill = colors_aty,
  main = "A. thalian",  # 隐藏标题 / Hide title
  category.names = names(aty_sets),
  # 优化标签位置，使其位于子集区域末端 / Optimize label positions to place them at the edge of the sets
  cat.pos = c(0,-72,-144,144,72),  
  alpha = venn_params$alpha,
  col = venn_params$col,
  cex = venn_params$cex,
  fontface = venn_params$fontface,
  cat.cex = venn_params$cat.cex,
  cat.fontface = venn_params$cat.fontface,
  cat.dist = venn_params$cat.dist,
  margin = venn_params$margin,
  sub.cex = venn_params$sub.cex,
  lwd = venn_params$lwd  # 应用更细的描边 / Apply thinner borders
)
grid.draw(venn_aty)
# 添加左上角的子图标签 'a' / Add subplot label 'a' in the upper left corner
grid.text("a", x = 0.05, y = 0.95, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()

# 第二个子图 (Ppa) / Second subplot (Ppa)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
# 针对4个集合的优化布局 / Optimized layout for 4 sets
venn_ppa <- venn.diagram(
  x = ppa_sets,
  filename = NULL,
  fill = colors_ppa,
  main = "P. patens",  # 隐藏标题 / Hide title
  category.names = names(ppa_sets),
  # 优化标签位置，使其位于子集区域末端 / Optimize label positions to place them at the edge of the sets
  cat.pos = c(-30, 30, -30, 30),  
  alpha = venn_params$alpha,
  col = venn_params$col,
  cex = venn_params$cex,
  fontface = venn_params$fontface,
  cat.cex = venn_params$cat.cex,
  cat.fontface = venn_params$cat.fontface,
  # cat.dist = venn_params$cat.dist,
  cat.dist = c(0.2,0.2,0.1,0.1),
  margin = venn_params$margin,
  sub.cex = venn_params$sub.cex,
  lwd = venn_params$lwd  # 应用更细的描边 / Apply thinner borders
)
grid.draw(venn_ppa)
# 添加左上角的子图标签 'b' / Add subplot label 'b' in the upper left corner
grid.text("b", x = 0.05, y = 0.95, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()

# 关闭PDF设备 / Close PDF device
dev.off()

cat(paste("\n\n改进后的维恩图已保存至", output_file, "\nImproved Venn diagram saved to", output_file))