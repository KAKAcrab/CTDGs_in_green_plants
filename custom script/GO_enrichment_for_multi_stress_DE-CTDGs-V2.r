#!/usr/bin/env Rscript
#########################################################################################################################
# 脚本名称: integrated_de_ctdg_go_enrichment.R
# 功能描述: 对多物种不同逆境条件下的DE-CTDGs进行统一的GO富集分析
# 优化特点: 1.过滤Significant值为0的条目 2.使用weight值替代KS值绘制点图 3.处理weight范围值 4.转换为对数值优化显示
#          5.修复classic和KS列显示，确保显示p-value值 6.保持结果与原始脚本完全一致
# 作者: YYH
# 日期: 2025-05-08
#########################################################################################################################

# 加载必要的R包 (Loading required packages)
suppressMessages({
  library(AnnotationHub)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)        # 用于数据转换和重塑
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
  library(graph)
  library(GO.db)
  library(SparseM)
  library(DOSE)
})

# 设置日志输出函数 (Setting up logging function)
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("开始运行改进版的DE-CTDGs富集分析脚本 v8.6 (Improved script started)")

# 配置并行计算 (Configure parallel computing)
num_cores <- detectCores() - 1  # 保留1个核心用于系统运行
if(num_cores < 1) num_cores <- 1 # 确保至少使用一个核心
log_message(paste0("使用 ", num_cores, " 个核心进行并行计算 (Using ", num_cores, " cores for parallel computing)"))
plan(multicore, workers = num_cores)

# 设置随机数种子以确保结果一致性 (Set random seed to ensure result consistency)
set.seed(12345)

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
    
    # 读取文件内容 (Read file contents)
    gene2go_data <- readMappings(file_path, sep = '\t')
    
    log_message(paste0("成功加载 ", length(gene2go_data), " 个基因的GO注释 (Successfully loaded GO annotations for ", 
                       length(gene2go_data), " genes)"))
    return(gene2go_data)
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

# 对指定条件的基因进行GO富集分析 - 保持原逻辑但添加p值提取
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
    # 使用统一的nodeSize参数以确保结果一致性
    GOdata <- new("topGOdata", ontology = ontology, allGenes = genelist,
                  annot = annFUN.gene2GO, gene2GO = gene2go,
                  nodeSize = 5) # 添加统一的节点大小阈值
    
    # 使用统一的统计检验方法 (Use unified statistical test methods)
    
    # 权重算法+Fisher检验 (weight algorithm + Fisher test)
    test.stat.weight <- new("weightCount", testStatistic = GOFisherTest, name="Fisher test", sigRatio = "ratio")
    resultWeight <- getSigGroups(GOdata, test.stat.weight)
    
    # 经典KS检验 (Classic KS test)
    test.stat.ks <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
    resultKS <- getSigGroups(GOdata, test.stat.ks)
    
    # 运行elim+ks算法 (Run elim+ks algorithm)
    elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    # 生成结果表格 (Generate result table)
    # 与原始代码保持一致的调用
    allRes <- GenTable(GOdata, classic=elim.ks, KS=resultKS, weight = resultWeight,
                       orderBy = "weight", ranksOf = "classic", topNodes = 15)
    
    # 确保结果数据框有足够的行 (Ensure the result data frame has enough rows)
    if(nrow(allRes) == 0) {
      log_message(paste0("警告: 条件 ", condition_name, " 的 ", ontology, " 富集分析没有找到显著的GO条目 (Warning: No significant GO terms found in ", 
                         ontology, " enrichment analysis for condition ", condition_name, ")"))
      return(data.frame())
    }
    
    # 添加条件信息和本体类型 (Add condition information and ontology type)
    allRes$Group <- condition_name
    allRes$ontology <- ontology
    
    # 将数值列转换为数值类型 (Convert numeric columns to numeric type)
    numeric_cols <- c("Annotated", "Significant", "Expected", "classic", "KS", "weight")
    for(col in numeric_cols) {
      if(col %in% colnames(allRes)) {
        allRes[[col]] <- as.numeric(as.character(allRes[[col]]))
      }
    }
    
    # 新增：获取p值并替换原始统计量
    # New: Get p-values and replace original statistics
    
    # 1. 获取GO ID列表
    go_ids <- allRes$GO.ID
    
    # 2. 运行额外的测试以获取p值
    classic_fisher_pvals <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    classic_ks_pvals <- runTest(GOdata, algorithm = "classic", statistic = "ks")
    
    # 3. 提取匹配的p值
    classic_fisher_scores <- score(classic_fisher_pvals)[go_ids]
    classic_ks_scores <- score(classic_ks_pvals)[go_ids]
    
    # 4. 替换到allRes数据框中
    allRes$classic_pvalue <- classic_fisher_scores
    allRes$KS_pvalue <- classic_ks_scores
    
    # 5. 如果需要，替换原列（注释掉这部分以保持原结果完全一致）
    allRes$classic <- classic_fisher_scores  # 此行将替换classic列为p值
    allRes$KS <- classic_ks_scores  # 此行将替换KS列为p值
    
    # 将结果转换为数据框 (Convert results to data frame)
    return(as.data.frame(allRes))
  }, error = function(e) {
    log_message(paste0("错误: 在进行条件 ", condition_name, " 的 ", ontology, " 富集分析时发生错误: ", e$message))
    return(data.frame())  # 返回空数据框
  })
}

# 并行执行所有条件的GO富集分析 (Perform GO enrichment analysis for all conditions in parallel)
perform_enrichment_for_all_conditions <- function(condition_genes, all_genes, gene2go, ontologies = c("BP", "MF", "CC")) {
  log_message(paste0("开始为所有条件执行富集分析，选择的本体类型: ", paste(ontologies, collapse = ", "), 
                     " (Starting enrichment analysis for all conditions with ontology types: ", 
                     paste(ontologies, collapse = ", "), ")"))
  
  all_results <- data.frame()
  
  # 对每种本体类型依次处理 (Process each ontology type sequentially)
  for(ontology in ontologies) {
    log_message(paste0("开始处理本体类型: ", ontology, " (Starting processing ontology type: ", ontology, ")"))
    
    # 使用future.apply来并行处理不同条件 (Use future.apply for parallel processing different conditions)
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
      all_results <- rbind(all_results, combined_results)
    } else {
      log_message(paste0("警告: 没有找到任何显著富集的 ", ontology, " GO条目 (Warning: No significant ", 
                         ontology, " GO terms found)"))
    }
  }
  
  if(nrow(all_results) > 0) {
    return(all_results)
  } else {
    log_message("警告: 所有本体类型均未发现显著富集的GO条目 (Warning: No significant GO terms found for any ontology type)")
    # 返回一个空的数据框但包含所有需要的列 (Return an empty data frame with all required columns)
    result_df <- data.frame(
      GO.ID = character(),
      Term = character(),
      Annotated = integer(),
      Significant = integer(),
      Expected = numeric(),
      classic = numeric(),
      KS = numeric(),
      weight = numeric(),
      Group = character(),
      ontology = character(),
      stringsAsFactors = FALSE
    )
    return(result_df)
  }
}

# 处理weight值中的范围值并计算对数转换
# Process range values in weight column and calculate logarithmic transformation
process_weight_values <- function(data) {
  # 确保weight列存在 (Ensure weight column exists)
  if(!"weight" %in% colnames(data)) {
    return(data)  # 如果不存在weight列，直接返回原数据
  }
  
  # 将weight值转换为字符串以便处理 (Convert weight values to string for processing)
  data$weight <- as.character(data$weight)
  
  # 处理含有"<"符号的范围值 (Process range values with "<" symbol)
  data$weight <- sapply(data$weight, function(w) {
    if(grepl("<", w)) {
      return(gsub("<", "", w))
    } else {
      return(w)
    }
  })
  
  # 将处理后的weight值转换为数值类型 (Convert processed weight values to numeric)
  data$weight <- as.numeric(data$weight)
  
  # 如果有NA值，替换为一个非常小的默认值 (If NA values exist, replace with a very small default value)
  data$weight[is.na(data$weight)] <- 1e-30
  
  # 计算-log10(weight)值用于点图 (Calculate -log10(weight) for plotting)
  data$weight_log <- -log10(data$weight)
  
  # 对classic和KS列也进行相同的处理 (Process classic and KS columns in the same way)
  if("classic" %in% colnames(data)) {
    data$classic <- as.character(data$classic)
    data$classic <- sapply(data$classic, function(c) {
      if(grepl("<", c)) {
        return(gsub("<", "", c))
      } else {
        return(c)
      }
    })
    data$classic <- as.numeric(data$classic)
    data$classic[is.na(data$classic)] <- 1e-30
    data$classic_log <- -log10(data$classic)
  }
  
  if("KS" %in% colnames(data)) {
    data$KS <- as.character(data$KS)
    data$KS <- sapply(data$KS, function(k) {
      if(grepl("<", k)) {
        return(gsub("<", "", k))
      } else {
        return(k)
      }
    })
    data$KS <- as.numeric(data$KS)
    data$KS[is.na(data$KS)] <- 1e-30
    data$KS_log <- -log10(data$KS)
  }
  
  return(data)
}

# 创建GO富集点图 - 为每种本体类型创建单独的点图
# (Create GO enrichment dot plot - create separate plots for each ontology type)
create_go_dotplot <- function(enrichment_data, output_dir, ontology = NULL, top_n = 15) {
  # 如果指定了本体类型，则筛选数据 (Filter data if ontology type is specified)
  if(!is.null(ontology)) {
    enrichment_data <- enrichment_data[enrichment_data$ontology == ontology, ]
    plot_title <- paste0(ontology, " Ontology Enrichment")
    output_file <- file.path(output_dir, paste0("DECTDGs_enrichment_of_", ontology, ".pdf"))
  } else {
    plot_title <- "GO Enrichment Analysis"
    output_file <- file.path(output_dir, "DECTDGs_enrichment_all_ontologies.pdf")
  }
  
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
  
  # 确保数值列是数值类型 (Ensure numeric columns are numeric type)
  numeric_cols <- c("Annotated", "Significant", "Expected", "classic", "KS", "weight")
  for(col in numeric_cols) {
    if(col %in% colnames(enrichment_data)) {
      enrichment_data[[col]] <- as.numeric(as.character(enrichment_data[[col]]))
    }
  }
  
  # 期望的条件顺序 (Expected condition order) - 根据实际情况可能需要修改
  # 自动检测所有条件并排序
  all_conditions <- unique(enrichment_data$Group)
  desired_order <- sort(all_conditions)  # 按字母顺序排序
  
  # 为每个条件选择前N个条目 (Select top N terms for each condition)
  conditions <- unique(enrichment_data$Group)
  # 确保条件按照期望顺序排序 (Ensure conditions are ordered as expected)
  conditions <- conditions[order(match(conditions, desired_order))]
  
  # 创建一个用于存储选中数据的数据框
  top_results <- data.frame()
  
  # 为每个条件处理数据
  for(condition in conditions) {
    # 提取当前条件的数据
    condition_data <- enrichment_data[enrichment_data$Group == condition, ]
    
    if(nrow(condition_data) > 0) {
      # 首先过滤掉Significant值为0的条目 (Filter out terms with Significant value of 0)
      condition_data <- condition_data %>% filter(Significant > 0)
      
      if(nrow(condition_data) > 0) {
        # 处理weight、classic和KS值中的范围值
        condition_data <- process_weight_values(condition_data)
        
        # 按weight值排序并选择前N个条目
        condition_top <- condition_data[order(condition_data$weight), ]
        condition_top <- head(condition_top, top_n)
        
        # 添加到结果数据框
        top_results <- rbind(top_results, condition_top)
      } else {
        log_message(paste0("警告: 条件 ", condition, " 过滤Significant=0后没有数据 (Warning: No data for condition ", 
                           condition, " after filtering Significant=0)"))
      }
    }
  }
  
  # 检查是否还有数据 (Check if there is still data)
  if(nrow(top_results) == 0) {
    log_message("警告: 筛选后没有数据，创建空白点图 (Warning: No data after filtering, creating empty dot plot)")
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = "No significant GO terms found") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_file, plot = p, width = 10, height = 10, units = "in")
    return(p)
  }
  
  # 创建点图所需的数据副本，以便修改显示格式而不影响原始数据
  plot_data <- top_results
  
  # 处理weight值并计算对数转换
  plot_data <- process_weight_values(plot_data)
  
  # ========= 处理重名GO条目的问题 =========
  # 仅针对显示标签，不修改原始数据
  # 使用GO.ID确保每个标签的唯一性
  
  # 1. 初始化绘图专用的显示标签列
  plot_data$DisplayTerm <- plot_data$Term
  
  # 2. 检查每个条件下是否有重复的Term
  for(condition in unique(plot_data$Group)) {
    condition_rows <- which(plot_data$Group == condition)
    if(length(condition_rows) > 0) {
      # 获取该条件下的所有Term
      condition_terms <- plot_data$Term[condition_rows]
      
      # 检查是否有重复的Term
      if(any(duplicated(condition_terms))) {
        # 处理每个重复的Term
        for(term in unique(condition_terms[duplicated(condition_terms)])) {
          # 找到有重复Term的所有行
          term_rows <- condition_rows[plot_data$Term[condition_rows] == term]
          
          # 给这些行的DisplayTerm添加GO.ID后缀
          for(row in term_rows) {
            # 使用完整GO ID作为标识
            full_id <- plot_data$GO.ID[row]
            # 更新显示标签
            plot_data$DisplayTerm[row] <- paste0(plot_data$Term[row], " [", full_id, "]")
          }
        }
      }
    }
  }
  # ==================
  
  # 将条件转换为因子，确保按期望顺序显示 
  # (Convert Group to factor to ensure display order as expected)
  plot_data$Group <- factor(plot_data$Group, levels = desired_order)
  
  # 将DisplayTerm转换为因子，确保按本体类型的顺序显示
  plot_data <- plot_data %>%
    arrange(ontology) %>%
    mutate(DisplayTerm = factor(DisplayTerm, levels = unique(DisplayTerm)))
  
  # 创建点图 - 使用weight_log替代KS (Create dot plot - use weight_log instead of KS)
  p <- ggplot(plot_data, aes(Group, DisplayTerm)) +
    geom_point(aes(color = weight_log, size = Significant)) +
    scale_color_gradient(
      name = "-log10(weight)")+
    scale_size_continuous(
      name = "Gene number"
    ) +
    theme_bw() +
    scale_color_gradient(low = '#6699CC', high = '#CC3333') + # 使用统一的颜色方案
    labs(x = NULL, y = NULL, title = plot_title) +
    guides(size = guide_legend(order = 1)) +
    theme(
      legend.direction = "vertical",
      legend.position = "right",
      axis.text.x = element_text(size = 8, colour = "black", face = "bold",angle = 0, hjust = .5),
      axis.text.y = element_text(size = 8, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    scale_y_discrete(position = "left")
  
  # 保存点图 (Save dot plot)
  ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in", dpi = 300)
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
  
  # 确保数值列是数值类型 (Ensure numeric columns are numeric type)
  numeric_cols <- c("Annotated", "Significant", "Expected", "classic", "KS", "weight")
  for(col in numeric_cols) {
    if(col %in% colnames(enrichment_data)) {
      enrichment_data[[col]] <- as.numeric(as.character(enrichment_data[[col]]))
    }
  }
  
  # 创建工作簿 (Create workbook)
  wb <- createWorkbook()
  
  # 添加汇总工作表 (Add summary worksheet)
  addWorksheet(wb, "Summary")
  
  # 首先过滤掉Significant值为0的记录
  filtered_enrichment_data <- enrichment_data %>% filter(Significant > 0)
  
  # 计算每个条件下各本体类型的GO条目数量 
  # (Calculate number of GO terms for each condition and ontology type)
  if(nrow(filtered_enrichment_data) > 0) {
    summary_data <- filtered_enrichment_data %>%
      group_by(Group, ontology) %>%
      summarise(Term_Count = n(), .groups = 'drop') %>%
      pivot_wider(names_from = ontology, values_from = Term_Count, values_fill = 0)
    
    # 写入汇总数据 (Write summary data)
    writeData(wb, "Summary", summary_data)
  } else {
    # 写入提示信息 (Write information)
    writeData(wb, "Summary", data.frame(Message = "No significant GO terms after filtering Significant=0"))
  }
  
  # 为每种本体类型创建单独的工作表 (Create separate worksheets for each ontology type)
  ontologies <- unique(enrichment_data$ontology)
  
  for(ontology in ontologies) {
    # 添加本体类型工作表 (Add ontology type worksheet)
    sheet_name <- paste0(ontology, "_Results")
    addWorksheet(wb, sheet_name)
    
    # 筛选当前本体类型的数据 (Filter data for current ontology type)
    ontology_data <- enrichment_data[enrichment_data$ontology == ontology, ]
    
    # 过滤Significant=0的记录
    ontology_data <- ontology_data %>% filter(Significant > 0)
    
    if(nrow(ontology_data) > 0) {
      # 处理weight值中的范围值
      ontology_data <- process_weight_values(ontology_data)
      
      # 创建一个包含所有条件TOP15结果的数据框 
      # (Create a data frame for top 15 results from all conditions)
      conditions <- unique(ontology_data$Group)
      # 自动排序条件 (Auto-sort conditions)
      conditions <- sort(conditions)
      
      ontology_top_results <- data.frame()
      
      for(condition in conditions) {
        # 提取条件数据并按weight值排序 (Extract condition data and sort by weight value)
        condition_data <- ontology_data[ontology_data$Group == condition, ]
        
        if(nrow(condition_data) > 0) {
          # 确保weight是数值类型
          condition_data$weight <- as.numeric(as.character(condition_data$weight))
          condition_data <- condition_data[order(condition_data$weight), ]
          
          # 选择前15个条目 (Select top 15 entries)
          if(nrow(condition_data) > 15) {
            condition_data <- head(condition_data, 15)
          }
          
          # 添加到结果数据框 (Add to results data frame)
          ontology_top_results <- rbind(ontology_top_results, condition_data)
        }
      }
      
      if(nrow(ontology_top_results) > 0) {
        # 按条件排序 (Sort by condition)
        ontology_top_results$Group <- factor(ontology_top_results$Group, levels = sort(unique(ontology_top_results$Group)))
        ontology_top_results$weight <- as.numeric(as.character(ontology_top_results$weight))
        ontology_top_results <- ontology_top_results[order(ontology_top_results$Group, ontology_top_results$weight), ]
        
        # 移除额外生成的对数列和p值列
        # Remove additional log columns and p-value columns
        extra_cols <- c("weight_log", "classic_log", "KS_log", "classic_pvalue", "KS_pvalue")
        ontology_top_results <- ontology_top_results[, !names(ontology_top_results) %in% extra_cols]
        
        # 写入数据 (Write data)
        writeData(wb, sheet_name, ontology_top_results)
        
        # 自动调整列宽 (Auto-adjust column widths)
        setColWidths(wb, sheet_name, cols = 1:ncol(ontology_top_results), widths = "auto")
      } else {
        # 没有数据写入提示信息 (No data, write information)
        writeData(wb, sheet_name, data.frame(Message = paste0("No significant ", ontology, " terms after filtering Significant=0")))
      }
    } else {
      # 没有数据写入提示信息 (No data, write information)
      writeData(wb, sheet_name, data.frame(Message = paste0("No significant ", ontology, " terms after filtering Significant=0")))
    }
  }
  
  # 添加包含所有结果的工作表 (Add worksheet with all results)
  addWorksheet(wb, "All_Results")
  
  # 过滤Significant=0的记录
  filtered_enrichment_data <- enrichment_data %>% filter(Significant > 0)
  
  if(nrow(filtered_enrichment_data) > 0) {
    # 处理weight值中的范围值
    filtered_enrichment_data <- process_weight_values(filtered_enrichment_data)
    
    # 移除额外生成的对数列和p值列
    # Remove additional log columns and p-value columns
    extra_cols <- c("weight_log", "classic_log", "KS_log", "classic_pvalue", "KS_pvalue")
    filtered_data_for_excel <- filtered_enrichment_data[, !names(filtered_enrichment_data) %in% extra_cols]
    
    # 按条件和本体类型排序 (Sort by condition and ontology type)
    filtered_data_for_excel$Group <- factor(filtered_data_for_excel$Group, levels = sort(unique(filtered_data_for_excel$Group)))
    filtered_data_for_excel$ontology <- factor(filtered_data_for_excel$ontology, levels = c("MF", "CC", "BP"))
    filtered_data_for_excel$weight <- as.numeric(as.character(filtered_data_for_excel$weight))
    filtered_data_for_excel <- filtered_data_for_excel[order(filtered_data_for_excel$Group, filtered_data_for_excel$ontology, filtered_data_for_excel$weight), ]
    
    # 写入所有数据 (Write all data)
    writeData(wb, "All_Results", filtered_data_for_excel)
    
    # 自动调整列宽 (Auto-adjust column widths)
    setColWidths(wb, "All_Results", cols = 1:ncol(filtered_data_for_excel), widths = "auto")
  } else {
    # 没有数据写入提示信息 (No data, write information)
    writeData(wb, "All_Results", data.frame(Message = "No significant GO terms after filtering Significant=0"))
  }
  
  # 保存工作簿 (Save workbook)
  saveWorkbook(wb, output_file, overwrite = TRUE)
  log_message(paste0("结果已保存到: ", output_file, " (Results saved to: ", output_file, ")"))
}

# 处理多物种 (Process multiple species)
process_multiple_species <- function(species_prefix, gene2go_file, all_genes_file, genes_file, expression_file, 
                                     output_dir, ontologies = c("BP", "MF", "CC")) {
  log_message(paste0("开始处理物种: ", species_prefix, " (Starting processing species: ", species_prefix, ")"))
  
  # 创建物种输出目录 (Create species output directory)
  species_output_dir <- file.path(output_dir, paste0(species_prefix, "_GO_enrichment"))
  dir.create(species_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 设置输出文件路径 (Set output file paths)
  output_xlsx <- file.path(species_output_dir, paste0(species_prefix, "_DECTDGs_GO_enrichment.xlsx"))
  
  # 读取输入文件 (Read input files)
  log_message("读取输入文件 (Reading input files)")
  
  # 读取基因到GO映射文件 (Read gene to GO mapping file)
  gene2go <- read_gene2go(gene2go_file)
  if(length(gene2go) == 0) {
    log_message("错误: 无法读取gene2go映射文件，跳过物种处理 (Error: Unable to read gene2go mapping file, skipping species processing)")
    return(FALSE)
  }
  
  # 读取所有基因列表 (Read all genes list)
  all_genes <- read_all_genes(all_genes_file)
  if(length(all_genes) == 0) {
    log_message("错误: 无法读取所有基因列表，跳过物种处理 (Error: Unable to read all genes list, skipping species processing)")
    return(FALSE)
  }
  
  # 读取基因簇成员数据 (Read gene cluster data)
  genes_data <- read_ctdg_genes(genes_file)
  if(nrow(genes_data) == 0) {
    log_message("错误: 无法读取基因簇成员数据，跳过物种处理 (Error: Unable to read gene cluster data, skipping species processing)")
    return(FALSE)
  }
  
  # 读取01分布矩阵 (Read binary expression matrix)
  binary_matrix <- read_expression_matrix(expression_file)
  if(nrow(binary_matrix) == 0) {
    log_message("错误: 无法读取01分布矩阵，跳过物种处理 (Error: Unable to read binary expression matrix, skipping species processing)")
    return(FALSE)
  }
  
  # 按条件分组DE-CTDGs (Group DE-CTDGs by condition)
  condition_genes <- group_dectdgs_by_condition(binary_matrix, genes_data)
  if(length(condition_genes) == 0) {
    log_message("错误: 没有找到任何条件下的DE-CTDGs，跳过物种处理 (Error: No DE-CTDGs found for any condition, skipping species processing)")
    return(FALSE)
  }
  
  # 执行GO富集分析，选择的本体类型 (Perform GO enrichment analysis for selected ontology types)
  enrichment_results <- perform_enrichment_for_all_conditions(condition_genes, all_genes, gene2go, ontologies)
  
  # 检查是否有结果 (Check if there are results)
  if(nrow(enrichment_results) > 0) {
    log_message("富集分析完成，开始创建可视化和保存结果 (Enrichment analysis completed, creating visualizations and saving results)")
    
    # 创建总体结果的点图 (Create overall result dot plot)
    overall_plot <- create_go_dotplot(enrichment_results, species_output_dir, NULL, 15)
    
    # 为每种本体类型创建单独的点图 (Create separate plots for each ontology type)
    for(ontology in ontologies) {
      ontology_plot <- create_go_dotplot(enrichment_results, species_output_dir, ontology, 15)
    }
    
    # 保存结果到Excel文件 (Save results to Excel file)
    save_to_excel(enrichment_results, output_xlsx)
    
    log_message(paste0("处理完成，结果已保存到: ", species_output_dir))
    return(TRUE)
  } else {
    log_message("警告: 没有找到任何显著富集的GO条目 (Warning: No significant GO terms found)")
    
    # 创建空的Excel文件 (Create empty Excel file)
    wb <- createWorkbook()
    addWorksheet(wb, "Info")
    writeData(wb, "Info", data.frame(Message = "No significant GO terms found"))
    saveWorkbook(wb, output_xlsx, overwrite = TRUE)
    
    # 创建空的点图 (Create empty plots)
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = "No significant GO terms found") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    output_file <- file.path(species_output_dir, "DECTDGs_enrichment_all_ontologies.pdf")
    ggsave(output_file, plot = p, width = 10, height = 10, units = "in")
    
    # 为每种本体类型创建空白点图
    for(ontology in ontologies) {
      output_file <- file.path(species_output_dir, paste0("DECTDGs_enrichment_of_", ontology, ".pdf"))
      ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in", dpi = 300)
    }
    
    log_message(paste0("处理完成，但没有找到显著结果，已创建空白输出文件: ", species_output_dir))
    return(FALSE)
  }
}

#########################################################################################################################
# 主函数 (Main function)
#########################################################################################################################

main <- function() {
  log_message("开始主函数执行 (Starting main function)")
  
  # 解析命令行参数 (Parse command line arguments)
  args <- commandArgs(trailingOnly = TRUE)
  
  # 默认参数 (Default parameters)
  species_prefix <- "Aty"  # 物种前缀，用于自动构建文件名
  gene2go_file <- paste0(species_prefix, "_gene2go")
  all_genes_file <- paste0(species_prefix, "_all_genes_list")
  genes_file <- paste0(species_prefix, "_processed_raw_genes_result")
  expression_file <- paste0(species_prefix, "_DE_CTDGs_matrix.csv")
  output_dir <- "GO_enrichment_results"
  # 默认分析所有本体类型 (Default to analyze all ontology types)
  # ontologies <- c("BP", "MF", "CC")
  ontologies <- c("BP")
  
  # 如果提供了命令行参数，则使用它们 (If command line arguments are provided, use them)
  if(length(args) >= 1) species_prefix <- args[1]
  if(length(args) >= 2) gene2go_file <- args[2]
  if(length(args) >= 3) all_genes_file <- args[3]
  if(length(args) >= 4) genes_file <- args[4]
  if(length(args) >= 5) expression_file <- args[5]
  if(length(args) >= 6) output_dir <- args[6]
  if(length(args) >= 7) {
    ont_arg <- args[7]
    # 解析本体类型参数 (Parse ontology type parameter)
    if(ont_arg != "") {
      # 支持多种分隔符号 (Support multiple separators)
      ontologies <- strsplit(ont_arg, "[,;|]")[[1]]
      ontologies <- toupper(trimws(ontologies)) # 转大写并去除空白 (Convert to uppercase and trim whitespace)
      # 验证本体类型 (Validate ontology types)
      valid_ontologies <- c("BP", "MF", "CC")
      ontologies <- ontologies[ontologies %in% valid_ontologies]
      if(length(ontologies) == 0) {
        log_message("警告: 未指定有效的本体类型，使用默认值 BP,MF,CC (Warning: No valid ontology types specified, using default BP,MF,CC)")
        ontologies <- c("BP", "MF", "CC")
      }
    }
  }
  
  # 创建总输出目录 (Create main output directory)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 检查是否处理多个物种 (Check if processing multiple species)
  if(length(args) >= 8 && args[8] == "multi") {
    # 处理多物种模式 (Multi-species mode)
    log_message("多物种模式已激活 (Multi-species mode activated)")
    
    # 解析物种列表 (Parse species list)
    species_list <- NULL
    if(length(args) >= 9) {
      species_list <- strsplit(args[9], "[,;|]")[[1]]
      species_list <- trimws(species_list)
    } else {
      # 尝试获取当前目录中所有符合条件的物种 (Try to get all species in current directory)
      matrix_files <- list.files(pattern = "_DE_CTDGs_matrix\\.csv$")
      species_list <- gsub("_DE_CTDGs_matrix\\.csv$", "", matrix_files)
    }
    
    if(length(species_list) == 0) {
      log_message("错误: 未找到任何物种数据，终止处理 (Error: No species data found, terminating)")
      return()
    }
    
    log_message(paste0("找到 ", length(species_list), " 个物种: ", paste(species_list, collapse = ", "), 
                       " (Found ", length(species_list), " species: ", paste(species_list, collapse = ", "), ")"))
    
    # 逐个处理每个物种 (Process each species)
    success_count <- 0
    for(species in species_list) {
      log_message(paste0("开始处理物种: ", species, " (Starting processing species: ", species, ")"))
      
      # 构建物种文件名 (Build species file names)
      sp_gene2go_file <- paste0(species, "_gene2go")
      sp_all_genes_file <- paste0(species, "_all_genes_list")
      sp_genes_file <- paste0(species, "_processed_raw_genes_result")
      sp_expression_file <- paste0(species, "_DE_CTDGs_matrix.csv")
      
      # 处理物种 (Process species)
      result <- process_multiple_species(species, sp_gene2go_file, sp_all_genes_file, 
                                         sp_genes_file, sp_expression_file, output_dir, ontologies)
      
      if(result) success_count <- success_count + 1
    }
    
    log_message(paste0("多物种处理完成，成功处理 ", success_count, " 个物种，共 ", length(species_list), " 个物种",
                       " (Multi-species processing completed, successfully processed ", success_count, 
                       " species out of ", length(species_list), " species)"))
  } else {
    # 单物种模式 (Single species mode)
    log_message(paste0("分析物种: ", species_prefix, 
                       "，使用本体类型: ", paste(ontologies, collapse = ", "), 
                       " (Analyzing species: ", species_prefix, 
                       ", using ontology types: ", paste(ontologies, collapse = ", "), ")"))
    
    # 处理单个物种 (Process single species)
    result <- process_multiple_species(species_prefix, gene2go_file, all_genes_file, 
                                       genes_file, expression_file, output_dir, ontologies)
    
    if(result) {
      log_message(paste0("物种 ", species_prefix, " 处理成功 (Species ", species_prefix, " processed successfully)"))
    } else {
      log_message(paste0("物种 ", species_prefix, " 处理完成，但可能没有发现显著富集 (Species ", 
                         species_prefix, " processing completed, but significant enrichment might not be found)"))
    }
  }
  
  log_message("脚本执行完成 (Script execution completed)")
}

# 执行主函数 (Execute main function)
main()