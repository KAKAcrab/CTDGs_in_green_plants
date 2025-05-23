#!/usr/bin/env Rscript
#########################################################################################################################
# 脚本名称: dectdg_condition_enrichment.R
# 功能描述: 基于不同逆境条件下的DE-CTDGs进行GO富集分析
# 日期: 2025-04-03
#########################################################################################################################

# 加载必要的R包 (Loading required packages)
suppressMessages({
  library(AnnotationHub)
  library(AnnotationDbi)
  library(dplyr)
  library(topGO)
  library(ggplot2)
  library(grid)
  library(Rgraphviz)
  library(gridExtra)
  library(openxlsx)
  library(parallel)
  library(doParallel)
  library(future)
  library(future.apply)
})

# 设置日志输出函数 (Setting up logging function)
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("开始运行逆境条件DE-CTDGs富集分析脚本 (Script started)")

# 配置并行计算 (Configure parallel computing)
num_cores <- detectCores() - 5  # 保留2个核心用于系统运行
if(num_cores < 1) num_cores <- 1 # 确保至少使用一个核心
log_message(paste0("使用 ", num_cores, " 个核心进行并行计算 (Using ", num_cores, " cores for parallel computing)"))
plan(multicore, workers = num_cores)

#########################################################################################################################
# 定义辅助函数 (Define helper functions)
#########################################################################################################################

# 读取基因到GO映射文件 (Read gene to GO mapping file)
read_gene2go <- function(file_path) {
  log_message(paste0("读取gene2go映射文件: ", file_path, " (Reading gene2go mapping file)"))
  
  tryCatch({
    # 检查文件是否存在
    if(!file.exists(file_path)) {
      stop(paste0("文件不存在: ", file_path))
    }
    
    # 读取文件内容
    gene2go_data <- read.table(file_path, sep = '\t', header = FALSE, 
                               stringsAsFactors = FALSE, quote = "", comment.char = "")
    colnames(gene2go_data) <- c("gene_id", "GO_terms")
    
    # 转换为topGO所需的映射格式
    gene2go <- list()
    for(i in 1:nrow(gene2go_data)) {
      gene_id <- gene2go_data$gene_id[i]
      go_terms <- unlist(strsplit(gene2go_data$GO_terms[i], ","))
      
      # 过滤掉空的GO项
      go_terms <- go_terms[go_terms != "" & !is.na(go_terms)]
      
      if(length(go_terms) > 0) {
        gene2go[[gene_id]] <- go_terms
      }
    }
    
    log_message(paste0("成功加载 ", length(gene2go), " 个基因的GO注释 (Successfully loaded GO annotations for ", length(gene2go), " genes)"))
    return(gene2go)
  }, error = function(e) {
    log_message(paste0("错误: 读取gene2go映射文件时出错: ", e$message))
    return(list())  # 返回空列表
  })
}

# 读取所有基因列表 (Read all genes list)
read_all_genes <- function(file_path) {
  log_message(paste0("读取所有基因列表: ", file_path, " (Reading all genes list)"))
  
  tryCatch({
    if(!file.exists(file_path)) {
      stop(paste0("文件不存在: ", file_path))
    }
    
    all_genes <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
    all_genes <- as.character(all_genes[,1])
    
    log_message(paste0("成功读取 ", length(all_genes), " 个基因 (Successfully read ", length(all_genes), " genes)"))
    return(all_genes)
  }, error = function(e) {
    log_message(paste0("错误: 读取所有基因列表时出错: ", e$message))
    return(character(0))  # 返回空字符向量
  })
}

# 读取基因簇成员数据 (Read gene cluster data)
read_ctdg_genes <- function(file_path) {
  log_message(paste0("读取CTDG基因数据: ", file_path, " (Reading CTDG gene data)"))
  
  tryCatch({
    if(!file.exists(file_path)) {
      stop(paste0("文件不存在: ", file_path))
    }
    
    # 检查第一行是否为标题行
    first_line <- readLines(file_path, n = 1)
    if(grepl("gene", tolower(first_line))) {
      skip_first_row <- TRUE
      log_message("检测到标题行，将跳过第一行 (Header detected, will skip first row)")
    } else {
      skip_first_row <- FALSE
    }
    
    # 读取文件
    if(skip_first_row) {
      genes_data <- read.table(file_path, sep = '\t', header = FALSE, skip = 1)
    } else {
      genes_data <- read.table(file_path, sep = '\t', header = FALSE)
    }
    
    colnames(genes_data) <- c("gene_id", "species_id", "chromosome_id", "cluster_id", "position_in_cluster")
    
    log_message(paste0("成功读取 ", nrow(genes_data), " 条基因簇成员数据 (Successfully read ", nrow(genes_data), " gene cluster records)"))
    return(genes_data)
  }, error = function(e) {
    log_message(paste0("错误: 读取CTDG基因数据时出错: ", e$message))
    return(data.frame())  # 返回空数据框
  })
}

# 读取表达矩阵 (Read expression matrix)
read_expression_matrix <- function(file_path) {
  log_message(paste0("读取表达矩阵: ", file_path, " (Reading expression matrix)"))
  
  tryCatch({
    if(!file.exists(file_path)) {
      stop(paste0("文件不存在: ", file_path))
    }
    
    # 读取表达矩阵
    expr_matrix <- read.csv(file_path, check.names = FALSE)
    
    # 检查数据格式
    if(ncol(expr_matrix) < 2) {
      stop("表达矩阵格式不正确，至少需要2列")
    }
    
    log_message(paste0("成功读取表达矩阵，包含 ", nrow(expr_matrix), " 个基因簇和 ", 
                       ncol(expr_matrix)-1, " 个条件 (Successfully read expression matrix with ", 
                       nrow(expr_matrix), " clusters and ", ncol(expr_matrix)-1, " conditions)"))
    
    return(expr_matrix)
  }, error = function(e) {
    log_message(paste0("错误: 读取表达矩阵时出错: ", e$message))
    return(data.frame())  # 返回空数据框
  })
}

# 按条件分组DE-CTDGs (Group DE-CTDGs by condition)
group_dectdgs_by_condition <- function(expr_matrix, genes_data) {
  log_message("按条件分组DE-CTDGs (Grouping DE-CTDGs by condition)")
  
  tryCatch({
    # 检查表达矩阵是否为空
    if(nrow(expr_matrix) == 0 || ncol(expr_matrix) <= 1) {
      stop("表达矩阵为空或只有一列")
    }
    
    # 获取条件列名（除了第一列以外的所有列）
    first_col_name <- colnames(expr_matrix)[1]
    condition_cols <- colnames(expr_matrix)[!colnames(expr_matrix) %in% c(first_col_name, "size")]
    
    # 如果没有找到条件列，终止处理
    if(length(condition_cols) == 0) {
      stop("未在表达矩阵中找到任何条件列")
    }
    
    log_message(paste0("找到 ", length(condition_cols), " 个条件: ", paste(condition_cols, collapse = ", ")))
    
    # 以条件为键创建列表
    condition_genes <- list()
    
    # 对每个条件，提取DE-CTDGs以及它们的成员基因
    for(condition in condition_cols) {
      # 提取当前条件下差异表达的基因簇
      de_clusters <- expr_matrix[expr_matrix[[condition]] == 1, first_col_name]
      
      if(length(de_clusters) > 0) {
        # 提取这些基因簇的所有成员基因
        de_genes <- genes_data[genes_data$cluster_id %in% de_clusters, "gene_id"]
        de_genes <- unique(de_genes)  # 去重
        
        condition_genes[[condition]] <- de_genes
        log_message(paste0("条件 ", condition, " 下有 ", length(de_clusters), " 个DE-CTDGs，包含 ", 
                           length(de_genes), " 个唯一基因 (Condition ", condition, " has ", 
                           length(de_clusters), " DE-CTDGs with ", length(de_genes), " unique genes)"))
      } else {
        condition_genes[[condition]] <- character(0)
        log_message(paste0("条件 ", condition, " 下没有DE-CTDGs (No DE-CTDGs found for condition ", condition, ")"))
      }
    }
    
    return(condition_genes)
  }, error = function(e) {
    log_message(paste0("错误: 按条件分组DE-CTDGs时出错: ", e$message))
    return(list())  # 返回空列表
  })
}

# 对指定条件的基因进行GO富集分析 (Perform GO enrichment analysis for genes in a specific condition)
perform_go_enrichment <- function(condition_name, condition_genes, all_genes, gene2go, ontology = "BP") {
  log_message(paste0("对条件 ", condition_name, " 的基因进行 ", ontology, " 本体富集分析 (Performing ", 
                     ontology, " ontology enrichment analysis for condition ", condition_name, ")"))
  
  # 检查条件是否为空
  if(length(condition_genes) == 0) {
    log_message(paste0("警告: 条件 ", condition_name, " 下没有基因，跳过富集分析 (Warning: No genes found for condition ", 
                       condition_name, ", skipping enrichment analysis)"))
    return(data.frame())
  }
  
  tryCatch({
    # 创建基因列表因子 (Create gene list factor)
    genenames <- names(gene2go)
    genelist <- factor(as.integer(genenames %in% condition_genes))
    names(genelist) <- genenames
    
    # 确保基因列表有至少一个位于GO注释中的基因 (Ensure gene list has at least one gene in GO annotations)
    if(sum(genelist == 1) == 0) {
      log_message(paste0("警告: 条件 ", condition_name, " 的基因与GO注释不匹配，跳过富集分析 (Warning: No genes for condition ", 
                         condition_name, " match GO annotations, skipping enrichment analysis)"))
      return(data.frame())
    }
    
    # 创建topGOdata对象 (Create topGOdata object)
    GOdata <- new("topGOdata", ontology = ontology, allGenes = genelist,
                  annot = annFUN.gene2GO, gene2GO = gene2go)
    
    # 执行Fisher精确检验 (Perform Fisher's exact test)
    test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    resultFisher <- getSigGroups(GOdata, test.stat)
    
    # 执行权重算法 (Perform weight algorithm)
    test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
    resultWeight <- getSigGroups(GOdata, test.stat)
    
    # 生成结果表格 (Generate result table)
    allRes <- GenTable(GOdata, classicFisher = resultFisher, weight = resultWeight,
                       orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30)
    
    # 确保结果数据框有足够的行 (Ensure the result data frame has enough rows)
    if(nrow(allRes) == 0) {
      log_message(paste0("警告: 条件 ", condition_name, " 的 ", ontology, " 富集分析没有找到显著的GO条目 (Warning: No significant GO terms found in ", 
                         ontology, " enrichment analysis for condition ", condition_name, ")"))
      return(data.frame())
    }
    
    # 添加条件信息 (Add condition information)
    allRes$Condition <- condition_name
    
    # 将classicFisher和weight列转换为数值 (Convert classicFisher and weight columns to numeric)
    allRes$classicFisher <- as.numeric(allRes$classicFisher)
    allRes$weight <- as.numeric(allRes$weight)
    
    return(as.data.frame(allRes))
  }, error = function(e) {
    log_message(paste0("错误: 在进行条件 ", condition_name, " 的 ", ontology, " 富集分析时发生错误: ", e$message))
    return(data.frame())  # 返回空数据框
  })
}

# 并行执行所有条件的GO富集分析 (Perform GO enrichment analysis for all conditions in parallel)
perform_enrichment_for_all_conditions <- function(condition_genes, all_genes, gene2go, ontology = "BP") {
  log_message(paste0("开始为所有条件执行 ", ontology, " 富集分析 (Starting ", ontology, " enrichment analysis for all conditions)"))
  
  # 使用future.apply来并行处理 (Use future.apply for parallel processing)
  results <- future_lapply(names(condition_genes), function(condition_name) {
    genes <- condition_genes[[condition_name]]
    res <- perform_go_enrichment(condition_name, genes, all_genes, gene2go, ontology)
    return(res)
  })
  
  # 为结果列表命名 (Name the results list)
  names(results) <- names(condition_genes)
  
  # 过滤空结果 (Filter empty results)
  filtered_results <- results[sapply(results, function(x) nrow(x) > 0)]
  
  if(length(filtered_results) > 0) {
    # 合并所有结果 (Combine all results)
    combined_results <- do.call(rbind, filtered_results)
    return(combined_results)
  } else {
    log_message(paste0("警告: 没有找到任何显著富集的 ", ontology, " GO条目 (Warning: No significant ", 
                       ontology, " GO terms found)"))
    # 返回一个空的数据框但包含所有需要的列 (Return an empty data frame with all required columns)
    result_df <- data.frame(
      GO.ID = character(),
      Term = character(),
      Annotated = integer(),
      Significant = integer(),
      Expected = numeric(),
      classicFisher = numeric(),
      weight = numeric(),
      Condition = character(),
      stringsAsFactors = FALSE
    )
    return(result_df)
  }
}

# 创建GO富集点图 (Create GO enrichment dot plot)
create_go_dotplot <- function(enrichment_data, output_file, title = "BP Ontology Enrichment") {
  log_message(paste0("创建GO富集点图: ", output_file, " (Creating GO enrichment dot plot)"))
  
  # 检查数据是否为空 (Check if data is empty)
  if(nrow(enrichment_data) == 0) {
    log_message("警告: 没有数据，创建空白点图 (Warning: No data, creating empty dot plot)")
    # 创建一个空的点图 (Create an empty plot)
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = "No significant GO terms found") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_file, plot = p, width = 10, height = 10, units = "in")
    return(p)
  }
  
  # 为每个条件选择前15个条目 (Select top 15 terms for each condition)
  conditions <- unique(enrichment_data$Condition)
  top_results <- data.frame()
  
  for(condition in conditions) {
    condition_data <- enrichment_data[enrichment_data$Condition == condition, ]
    
    if(nrow(condition_data) > 0) {
      # 按p值排序并选择前15个条目 (Sort by p-value and select top 15 terms)
      condition_top <- condition_data[order(condition_data$classicFisher), ]
      condition_top <- head(condition_top, 15)
      top_results <- rbind(top_results, condition_top)
    }
  }
  
  # 为重复的Terms添加一个显示用的列，但不包含GO ID
  # (Add a display column for duplicate terms, without GO ID)
  top_results$DisplayTerm <- top_results$Term
  
  # 检查是否有重复的术语 (Check for duplicate terms)
  term_count <- table(top_results$Term)
  duplicate_terms <- names(term_count[term_count > 1])
  
  # 对于在不同条件下出现的相同术语，我们不做特殊处理
  # 因为它们会按条件分组显示，不会在图表中产生混淆
  # (For terms appearing across different conditions, no special handling needed
  # as they will be grouped by condition in the plot)
  
  # 将DisplayTerm转换为因子，确保它们以特定顺序显示 
  # (Convert DisplayTerm to factor to ensure specific display order)
  # 首先按条件分组，然后按p值排序
  ordered_terms <- c()
  for(condition in conditions) {
    condition_data <- top_results[top_results$Condition == condition, ]
    condition_data <- condition_data[order(condition_data$classicFisher), ]
    ordered_terms <- c(ordered_terms, condition_data$DisplayTerm)
  }
  unique_ordered_terms <- unique(ordered_terms)
  top_results$DisplayTerm <- factor(top_results$DisplayTerm, levels = rev(unique_ordered_terms))
  
  # 创建点图 (Create dot plot)
  p <- ggplot(top_results, aes(x = Condition, y = DisplayTerm)) +
    geom_point(aes(color = -log10(classicFisher), size = Significant)) +
    theme_bw() +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(p-value)") +
    scale_size(name = "Gene Count") +
    labs(x = NULL, y = NULL, title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey95")
    )
  
  # 保存点图 (Save dot plot)
  ggsave(output_file, plot = p, width = 10, height = 12, units = "in", dpi = 300)
  log_message(paste0("点图已保存到: ", output_file, " (Plot saved to: ", output_file, ")"))
  
  return(p)
}

# 保存结果到Excel文件 (Save results to Excel file)
save_to_excel <- function(enrichment_data, output_file) {
  log_message(paste0("保存结果到Excel文件: ", output_file, " (Saving results to Excel file)"))
  
  # 检查数据是否为空 (Check if data is empty)
  if(nrow(enrichment_data) == 0) {
    log_message("警告: 没有数据可以保存 (Warning: No data to save)")
    # 创建一个空的工作簿 (Create an empty workbook)
    wb <- createWorkbook()
    addWorksheet(wb, "Info")
    writeData(wb, "Info", data.frame(Message = "No significant GO terms found"))
    saveWorkbook(wb, output_file, overwrite = TRUE)
    return()
  }
  
  # 创建工作簿 (Create workbook)
  wb <- createWorkbook()
  
  # 添加汇总工作表 (Add summary worksheet)
  addWorksheet(wb, "Summary")
  
  # 计算每个条件的GO条目数量 (Calculate number of GO terms for each condition)
  summary_data <- table(enrichment_data$Condition)
  summary_df <- data.frame(
    Condition = names(summary_data),
    GO_Term_Count = as.numeric(summary_data)
  )
  
  # 写入汇总数据 (Write summary data)
  writeData(wb, "Summary", summary_df)
  
  # 创建一个包含所有条件TOP30结果的数据框 (Create a data frame for top 30 results from all conditions)
  conditions <- unique(enrichment_data$Condition)
  all_top_results <- data.frame()
  
  for(condition in conditions) {
    # 提取条件数据并按p值排序 (Extract condition data and sort by p-value)
    condition_data <- enrichment_data[enrichment_data$Condition == condition, ]
    condition_data <- condition_data[order(condition_data$classicFisher), ]
    
    # 选择前30个条目 (Select top 30 entries)
    if(nrow(condition_data) > 30) {
      condition_data <- head(condition_data, 30)
    }
    
    # 添加到结果数据框 (Add to results data frame)
    all_top_results <- rbind(all_top_results, condition_data)
  }
  
  # 按条件排序 (Sort by condition)
  all_top_results <- all_top_results[order(all_top_results$Condition, all_top_results$classicFisher), ]
  
  # 添加结果工作表 (Add results worksheet)
  addWorksheet(wb, "Results")
  
  # 写入数据 (Write data)
  writeData(wb, "Results", all_top_results)
  
  # 自动调整列宽 (Auto-adjust column widths)
  setColWidths(wb, "Results", cols = 1:ncol(all_top_results), widths = "auto")
  
  # 保存工作簿 (Save workbook)
  saveWorkbook(wb, output_file, overwrite = TRUE)
  log_message(paste0("结果已保存到: ", output_file, " (Results saved to: ", output_file, ")"))
}

#########################################################################################################################
# 主函数 (Main function)
#########################################################################################################################

main <- function() {
  log_message("开始主函数执行 (Starting main function)")
  
  # 解析命令行参数 (Parse command line arguments)
  args <- commandArgs(trailingOnly = TRUE)
  
  # 默认参数 (Default parameters)
  binary_matrix_file <- "Aty_DE_CTDGs_matrix.csv"
  genes_file <- "Aty_processed_raw_genes_result"
  gene2go_file <- "Aty_gene2go"
  all_genes_file <- "Aty_all_genes_list"
  output_dir <- "Aty_GO_enrichment_results"
  
  # 如果提供了命令行参数，则使用它们 (If command line arguments are provided, use them)
  if(length(args) >= 1) binary_matrix_file <- args[1]
  if(length(args) >= 2) genes_file <- args[2]
  if(length(args) >= 3) gene2go_file <- args[3]
  if(length(args) >= 4) all_genes_file <- args[4]
  if(length(args) >= 5) output_dir <- args[5]
  
  # 创建输出目录 (Create output directory)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 设置输出文件路径 (Set output file paths)
  output_plot <- file.path(output_dir, "DECTDGs_condition_BP_enrichment.pdf")
  output_xlsx <- file.path(output_dir, "DECTDGs_condition_BP_enrichment.xlsx")
  
  # 读取输入文件 (Read input files)
  log_message("读取输入文件 (Reading input files)")
  
  # 读取基因到GO映射文件 (Read gene to GO mapping file)
  gene2go <- read_gene2go(gene2go_file)
  if(length(gene2go) == 0) {
    log_message("错误: 无法读取gene2go映射文件，终止处理 (Error: Unable to read gene2go mapping file, terminating)")
    return()
  }
  
  # 读取所有基因列表 (Read all genes list)
  all_genes <- read_all_genes(all_genes_file)
  if(length(all_genes) == 0) {
    log_message("错误: 无法读取所有基因列表，终止处理 (Error: Unable to read all genes list, terminating)")
    return()
  }
  
  # 读取基因簇成员数据 (Read gene cluster data)
  genes_data <- read_ctdg_genes(genes_file)
  if(nrow(genes_data) == 0) {
    log_message("错误: 无法读取基因簇成员数据，终止处理 (Error: Unable to read gene cluster data, terminating)")
    return()
  }
  
  # 读取01分布矩阵 (Read binary expression matrix)
  binary_matrix <- read_expression_matrix(binary_matrix_file)
  if(nrow(binary_matrix) == 0) {
    log_message("错误: 无法读取01分布矩阵，终止处理 (Error: Unable to read binary expression matrix, terminating)")
    return()
  }
  
  # 按条件分组DE-CTDGs (Group DE-CTDGs by condition)
  condition_genes <- group_dectdgs_by_condition(binary_matrix, genes_data)
  if(length(condition_genes) == 0) {
    log_message("错误: 没有找到任何条件下的DE-CTDGs，终止处理 (Error: No DE-CTDGs found for any condition, terminating)")
    return()
  }
  
  # 执行GO富集分析，只关注BP本体 (Perform GO enrichment analysis, focusing on BP ontology)
  enrichment_results <- perform_enrichment_for_all_conditions(condition_genes, all_genes, gene2go, "BP")
  
  # 检查是否有结果 (Check if there are results)
  if(nrow(enrichment_results) > 0) {
    # 创建可视化 (Create visualization)
    bp_plot <- create_go_dotplot(enrichment_results, output_plot, "Biological Process Enrichment for DE-CTDGs")
    
    # 保存结果到Excel文件 (Save results to Excel file)
    save_to_excel(enrichment_results, output_xlsx)
    
    log_message(paste0("处理完成，结果已保存到:\n", 
                       output_plot, "\n", 
                       output_xlsx))
  } else {
    log_message("警告: 没有找到任何显著富集的GO条目 (Warning: No significant GO terms found)")
    
    # 创建空的点图 (Create empty plot)
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = "No significant GO terms found") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_plot, plot = p, width = 10, height = 10, units = "in")
    
    # 创建空的Excel文件 (Create empty Excel file)
    wb <- createWorkbook()
    addWorksheet(wb, "Info")
    writeData(wb, "Info", data.frame(Message = "No significant GO terms found"))
    saveWorkbook(wb, output_xlsx, overwrite = TRUE)
    
    log_message(paste0("处理完成，但没有找到显著结果，已创建空白输出文件:\n", 
                       output_plot, "\n", 
                       output_xlsx))
  }
  
  log_message("脚本执行完成 (Script execution completed)")
}

# 执行主函数 (Execute main function)
main()