#!/usr/bin/env Rscript
#########################################################################################################################
# 功能描述: 基于不同大小的CTDGs进行GO富集分析
# 优化特点: 1.过滤Significant值为0的条目 2.使用weight值替代KS值绘制点图 3.处理weight范围值 4.转换为对数值优化显示
#          5.修复classic和KS列显示，确保显示p-value值 6.保持结果与原始脚本完全一致
# 日期: 2025-05-08
#########################################################################################################################

# 加载必要的R包 (Loading required packages)
suppressMessages({
  library(AnnotationHub)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
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

log_message("开始运行CTDG大小分组富集分析脚本 (Script started)")

# 配置并行计算 (Configure parallel computing)
num_cores <- detectCores() - 3  # 保留3个核心用于系统运行
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

# 读取并处理CTDG基因数据 (Read and process CTDG gene data)
read_ctdg_genes <- function(file_path) {
  log_message(paste0("读取CTDG基因数据: ", file_path, " (Reading CTDG gene data)"))
  
  tryCatch({
    # 检查文件是否存在
    if(!file.exists(file_path)) {
      stop(paste0("文件不存在: ", file_path))
    }
    
    genes_data <- read.table(file_path, sep = '\t', header = FALSE)
    colnames(genes_data) <- c("gene_id", "species_id", "chromosome_id", "cluster_id", "position_in_cluster")
    
    log_message(paste0("成功读取 ", nrow(genes_data), " 条基因簇成员数据 (Successfully read ", nrow(genes_data), " gene cluster records)"))
    return(genes_data)
  }, error = function(e) {
    log_message(paste0("错误: 读取CTDG基因数据时出错: ", e$message))
    return(data.frame())  # 返回空数据框
  })
}

# 读取并处理CTDG簇信息 (Read and process CTDG cluster information)
read_ctdg_numbers <- function(file_path) {
  log_message(paste0("读取CTDG簇信息: ", file_path, " (Reading CTDG cluster information)"))
  
  tryCatch({
    # 检查文件是否存在
    if(!file.exists(file_path)) {
      stop(paste0("文件不存在: ", file_path))
    }
    
    numbers_data <- read.table(file_path, sep = '\t', header = FALSE)
    colnames(numbers_data) <- c("species_id", "chromosome_id", "cluster_id", "member_count")
    
    log_message(paste0("成功读取 ", nrow(numbers_data), " 条CTDG簇信息 (Successfully read ", nrow(numbers_data), " CTDG cluster records)"))
    return(numbers_data)
  }, error = function(e) {
    log_message(paste0("错误: 读取CTDG簇信息时出错: ", e$message))
    return(data.frame())  # 返回空数据框
  })
}

# 读取所有基因列表 (Read all genes list)
read_all_genes <- function(file_path) {
  log_message(paste0("读取所有基因列表: ", file_path, " (Reading all genes list)"))
  
  tryCatch({
    # 检查文件是否存在
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

# 根据簇大小对基因进行分组 (Group genes by cluster size)
group_genes_by_cluster_size <- function(genes_data, numbers_data) {
  log_message("根据基因簇大小对基因进行分组 (Grouping genes by cluster size)")
  
  tryCatch({
    # 合并基因数据和簇大小信息
    # Merge gene data with cluster size information
    gene_with_size <- genes_data %>%
      left_join(numbers_data, by = c("species_id", "chromosome_id", "cluster_id"))
    
    # 根据簇大小将基因分为四组
    # Group genes into four categories based on cluster size
    group1 <- gene_with_size %>% filter(member_count == 2) %>% pull(gene_id) %>% unique()
    group2 <- gene_with_size %>% filter(member_count >= 3 & member_count <= 5) %>% pull(gene_id) %>% unique()
    group3 <- gene_with_size %>% filter(member_count >= 6 & member_count <= 10) %>% pull(gene_id) %>% unique()
    group4 <- gene_with_size %>% filter(member_count > 10) %>% pull(gene_id) %>% unique()
    
    groups <- list(
      "size_2" = group1,
      "size_3_5" = group2,
      "size_6_10" = group3,
      "size_gt_10" = group4
    )
    
    # 输出每组基因的数量
    # Output the number of genes in each group
    for(i in 1:length(groups)) {
      log_message(paste0("组 ", names(groups)[i], " 包含 ", length(groups[[i]]), " 个基因 (Group ", 
                         names(groups)[i], " contains ", length(groups[[i]]), " genes)"))
    }
    
    return(groups)
  }, error = function(e) {
    log_message(paste0("错误: 基因簇大小分组时出错: ", e$message))
    return(list())  # 返回空列表
  })
}

# 对指定组的基因进行GO富集分析 - 保持原逻辑但添加p值提取
perform_go_enrichment <- function(group_name, group_genes, all_genes, gene2go, ontology = "BP") {
  log_message(paste0("对组 ", group_name, " 的基因进行 ", ontology, " 本体富集分析 (Performing ", 
                     ontology, " ontology enrichment analysis for group ", group_name, ")"))
  
  # 检查组是否为空
  if(length(group_genes) == 0) {
    log_message(paste0("警告: 组 ", group_name, " 内没有基因，跳过富集分析 (Warning: No genes in group ", 
                       group_name, ", skipping enrichment analysis)"))
    return(data.frame())
  }
  
  tryCatch({
    # 创建基因列表因子 (Create gene list factor)
    genenames <- names(gene2go)
    genelist <- factor(as.integer(genenames %in% group_genes))
    names(genelist) <- genenames
    
    # 确保基因列表有至少一个位于GO注释中的基因 (Ensure gene list has at least one gene in GO annotations)
    if(sum(genelist == 1) == 0) {
      log_message(paste0("警告: 组 ", group_name, " 的基因与GO注释不匹配，跳过富集分析 (Warning: No genes in group ", 
                         group_name, " match GO annotations, skipping enrichment analysis)"))
      return(data.frame())
    }
    
    # 检查因子级别 (Check factor levels)
    if(length(levels(genelist)) < 2) {
      log_message(paste0("警告: 因子级别不足，跳过 ", ontology, " 富集分析 (Warning: Insufficient factor levels, skipping ", 
                         ontology, " enrichment analysis)"))
      return(data.frame())
    }
    
    # 创建topGOdata对象 (Create topGOdata object)
    GOdata <- new("topGOdata", ontology = ontology, allGenes = genelist,
                  annot = annFUN.gene2GO, gene2GO = gene2go)
    
    # # 经典Fisher检验 (Classic Fisher test)
    # test.stat.fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
    # resultFisher <- getSigGroups(GOdata, test.stat.fisher)
    # 
    # # elim算法+KS检验 (elim algorithm + KS test)
    # test.stat.elim <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
    # resultElim <- getSigGroups(GOdata, test.stat.elim)
    
    # 权重算法+Fisher检验 (weight algorithm + Fisher test)
    test.stat.weight <- new("weightCount", testStatistic = GOFisherTest, name="Fisher test", sigRatio = "ratio")
    resultWeight <- getSigGroups(GOdata, test.stat.weight)
    
    # 经典KS检验 (Classic KS test)
    test.stat.ks <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
    resultKS <- getSigGroups(GOdata, test.stat.ks)
    
    # 运行elim+ks算法 (Run elim+ks algorithm)
    elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
    
    # 保留前15个条目，与原始脚本保持一致
    allRes <- GenTable(GOdata, classic=elim.ks, KS=resultKS, weight = resultWeight,
                       orderBy = "weight", ranksOf = "classic", topNodes = 15)
    
    # 确保结果数据框有足够的行 (Ensure the result data frame has enough rows)
    if(nrow(allRes) == 0) {
      log_message(paste0("警告: 组 ", group_name, " 的 ", ontology, " 富集分析没有找到显著的GO条目 (Warning: No significant GO terms found in ", 
                         ontology, " enrichment analysis for group ", group_name, ")"))
      return(data.frame())
    }
    
    # 添加组信息和本体类型 (Add group information and ontology type)
    allRes$Group <- group_name
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
    
    # 5. 替换原列 - 将classic和KS列的值替换为p值
    allRes$classic <- classic_fisher_scores
    allRes$KS <- classic_ks_scores
    
    # 将结果转换为数据框 (Convert results to data frame)
    return(as.data.frame(allRes))
  }, error = function(e) {
    log_message(paste0("错误: 在进行组 ", group_name, " 的 ", ontology, " 富集分析时发生错误: ", e$message))
    return(data.frame())  # 返回空数据框
  })
}

# 执行所有本体的GO富集分析 (Perform GO enrichment analysis for all ontologies)
perform_all_ontologies <- function(group_name, group_genes, all_genes, gene2go) {
  log_message(paste0("为组 ", group_name, " 执行所有本体的GO富集分析 (Performing GO enrichment for all ontologies for group ", group_name, ")"))
  
  # 检查组是否为空 (Check if group is empty)
  if(length(group_genes) == 0) {
    log_message(paste0("警告: 组 ", group_name, " 内没有基因，跳过所有富集分析 (Warning: No genes in group ", group_name, ", skipping all enrichment analyses)"))
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
  
  # 对三种不同的本体执行富集分析 (Perform enrichment for three different ontologies)
  ontologies <- c("MF", "BP", "CC")
  
  results <- list()
  
  for(ont in ontologies) {
    res <- perform_go_enrichment(group_name, group_genes, all_genes, gene2go, ont)
    if(nrow(res) > 0) {
      results[[ont]] <- res
    }
  }
  
  # 合并结果 (Combine results)
  if(length(results) > 0) {
    combined_results <- do.call(rbind, results)
    return(combined_results)
  } else {
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

# 并行执行所有组的GO富集分析 (Perform GO enrichment analysis for all groups in parallel)
perform_enrichment_for_all_groups <- function(gene_groups, all_genes, gene2go) {
  log_message("开始为所有组执行GO富集分析 (Starting GO enrichment analysis for all groups)")
  
  # 使用future.apply来并行处理 (Use future.apply for parallel processing)
  results <- future_lapply(names(gene_groups), function(group_name) {
    group_genes <- gene_groups[[group_name]]
    res <- perform_all_ontologies(group_name, group_genes, all_genes, gene2go)
    return(res)
  })
  
  # 合并所有结果 (Combine all results)
  filtered_results <- results[sapply(results, function(x) nrow(x) > 0)]
  
  if(length(filtered_results) > 0) {
    combined_results <- do.call(rbind, filtered_results)
    return(combined_results)
  } else {
    log_message("警告: 没有找到任何显著富集的GO条目 (Warning: No significant GO terms found)")
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

# 创建GO富集点图 (Create GO enrichment dot plot)
create_go_dotplot <- function(enrichment_data, ont, output_file) {
  log_message(paste0("为 ", ont, " 本体创建点图 (Creating dot plot for ", ont, " ontology)"))
  
  # 检查数据是否为空 (Check if data is empty)
  if(nrow(enrichment_data) == 0) {
    log_message(paste0("警告: 没有 ", ont, " 本体的数据，跳过创建点图 (Warning: No data for ", ont, " ontology, skipping dot plot creation)"))
    # 创建一个空的点图 (Create an empty plot)
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = paste0("No significant ", ont, " terms found")) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in")
    return(p)
  }
  
  # 筛选指定本体的数据 (Filter data for specified ontology)
  ont_data <- enrichment_data %>% 
    filter(ontology == ont)
  
  # 确保数值列是数值类型 (Ensure numeric columns are numeric type)
  numeric_cols <- c("Annotated", "Significant", "Expected", "classic", "KS", "weight")
  for(col in numeric_cols) {
    if(col %in% colnames(ont_data)) {
      ont_data[[col]] <- as.numeric(as.character(ont_data[[col]]))
    }
  }
  
  # 如果筛选后数据为空，则返回一个空的点图 (If filtered data is empty, return an empty plot)
  if(nrow(ont_data) == 0) {
    log_message(paste0("警告: 筛选后没有 ", ont, " 本体的数据，跳过创建点图 (Warning: No filtered data for ", ont, " ontology, skipping dot plot creation)"))
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = paste0("No significant ", ont, " terms found")) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in")
    return(p)
  }
  
  # 为每个组选择前15个条目 (Select top 15 terms for each group)
  groups <- unique(ont_data$Group)
  top_results <- data.frame()
  
  for(group in groups) {
    group_data <- ont_data %>% filter(Group == group)
    
    if(nrow(group_data) > 0) {
      # 首先过滤掉Significant值为0的条目 (Filter out terms with Significant value of 0)
      group_data <- group_data %>% filter(Significant > 0)
      
      if(nrow(group_data) > 0) {
        # 处理weight值中的范围值
        # 将所有weight值转换为字符串以便处理
        group_data$weight <- as.character(group_data$weight)
        
        # 处理含有"<"符号的范围值
        group_data$weight <- sapply(group_data$weight, function(w) {
          if(grepl("<", w)) {
            return(gsub("<", "", w))
          } else {
            return(w)
          }
        })
        
        # 将处理后的weight值转换为数值类型
        group_data$weight <- as.numeric(group_data$weight)
        
        # 如果有NA值，替换为一个非常小的默认值
        group_data$weight[is.na(group_data$weight)] <- 1e-30
        
        # 按weight值排序并选择前15个条目
        group_top <- group_data %>% 
          arrange(weight) %>% 
          head(15)
        
        # 添加到结果数据框
        top_results <- rbind(top_results, group_top)
      } else {
        log_message(paste0("警告: 组 ", group, " 过滤Significant=0后没有数据 (Warning: No data for group ", 
                           group, " after filtering Significant=0)"))
      }
    }
  }
  
  # 如果没有结果，则返回空的点图 (If no results, return empty plot)
  if(nrow(top_results) == 0) {
    log_message(paste0("警告: 没有选出任何 ", ont, " 本体的条目，跳过创建点图 (Warning: No terms selected for ", ont, " ontology, skipping dot plot creation)"))
    p <- ggplot() + 
      theme_void() + 
      annotate("text", x = 0, y = 0, label = paste0("No significant ", ont, " terms found")) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in")
    return(p)
  }
  
  # 处理weight值并计算对数转换
  top_results <- process_weight_values(top_results)
  
  # 将Term转换为因子以保持顺序 (Convert Term to factor to maintain order)
  top_results$Term <- factor(top_results$Term, levels = unique(top_results$Term))
  
  # 确保Group列的顺序 (Ensure proper order of Group column)
  desired_order <- c("size_2", "size_3_5", "size_6_10", "size_gt_10")
  valid_groups <- desired_order[desired_order %in% unique(top_results$Group)]
  top_results$Group <- factor(top_results$Group, levels = valid_groups)
  
  # 创建点图 - 使用weight的对数值 (Create dot plot - use log of weight value for color)
  p <- ggplot(top_results, aes(Group, Term)) +
    geom_point(aes(color = weight_log, size = Significant)) +
    theme_bw() +
    scale_size_continuous(name = "Gene number")+
    scale_color_gradient(low = '#6699CC', high = '#CC3333', name = "-log10(weight)") +
    labs(x = NULL, y = NULL, title = paste0(ont, " Ontology Enrichment")) +
    guides(size = guide_legend(order = 1)) +
    theme(
      legend.direction = "vertical",
      legend.position = "right",
      axis.text.x = element_text(size = 8, colour = "black", face = "bold"),
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
    log_message("警告: 没有数据可以保存到Excel文件 (Warning: No data to save to Excel file)")
    # 创建一个包含提示信息的工作簿 (Create a workbook with information)
    wb <- createWorkbook()
    addWorksheet(wb, "Info")
    writeData(wb, "Info", data.frame(Message = "No significant GO terms found"))
    saveWorkbook(wb, output_file, overwrite = TRUE)
    return(invisible(NULL))
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
  
  # 为每个本体创建工作表 (Create worksheet for each ontology)
  ontologies <- c("MF", "BP", "CC")
  
  for(ont in ontologies) {
    # 筛选本体数据 (Filter data for ontology)
    ont_data <- enrichment_data %>% 
      filter(ontology == ont)
    
    if(nrow(ont_data) > 0) {
      # 添加工作表 (Add worksheet)
      addWorksheet(wb, ont)
      
      # 为每个组选择前15个条目 (Select top 15 terms for each group)
      groups <- unique(ont_data$Group)
      top_results <- data.frame()
      
      for(group in groups) {
        group_data <- ont_data %>% filter(Group == group)
        
        # 筛选出Significant值大于0的条目 (Filter terms with Significant value greater than 0)
        group_data <- group_data %>% filter(Significant > 0)
        
        # 如果过滤后还有数据，则处理weight值和排序 (If there is still data after filtering, process weight values and sort)
        if(nrow(group_data) > 0) {
          # 处理weight值中的范围值
          # 将所有weight值转换为字符串以便处理
          group_data$weight <- as.character(group_data$weight)
          
          # 处理含有"<"符号的范围值
          group_data$weight <- sapply(group_data$weight, function(w) {
            if(grepl("<", w)) {
              return(gsub("<", "", w))
            } else {
              return(w)
            }
          })
          
          # 将处理后的weight值转换为数值类型
          group_data$weight <- as.numeric(group_data$weight)
          
          # 如果有NA值，替换为一个非常小的默认值
          group_data$weight[is.na(group_data$weight)] <- 1e-30
          
          # 按weight值排序并选择前15个条目
          group_top <- group_data %>% 
            arrange(weight) %>% 
            head(15)
          
          # 移除额外生成的对数列和p值列
          extra_cols <- c("weight_log", "classic_log", "KS_log", "classic_pvalue", "KS_pvalue") 
          if(any(extra_cols %in% colnames(group_top))) {
            group_top <- group_top[, !names(group_top) %in% extra_cols]
          }
          
          top_results <- rbind(top_results, group_top)
        } else {
          log_message(paste0("警告: 组 ", group, " 过滤Significant=0后没有数据 (Warning: No data for group ", 
                             group, " after filtering Significant=0)"))
        }
      }
      
      # 确保Group列的顺序 (Ensure proper order of Group column)
      if(nrow(top_results) > 0) {
        desired_order <- c("size_2", "size_3_5", "size_6_10", "size_gt_10")
        valid_groups <- desired_order[desired_order %in% unique(top_results$Group)]
        top_results$Group <- factor(top_results$Group, levels = valid_groups)
        top_results$weight <- as.numeric(as.character(top_results$weight))
        top_results <- top_results[order(top_results$Group, top_results$weight), ]
        
        # 写入数据 (Write data)
        writeData(wb, ont, top_results)
        
        # 设置列宽 (Set column widths)
        setColWidths(wb, ont, cols = 1:ncol(top_results), widths = "auto")
      } else {
        # 添加空的工作表 (Add empty worksheet)
        writeData(wb, ont, data.frame(Message = paste0("No significant ", ont, " terms found")))
      }
    } else {
      # 添加空的工作表 (Add empty worksheet)
      addWorksheet(wb, ont)
      writeData(wb, ont, data.frame(Message = paste0("No significant ", ont, " terms found")))
    }
  }
  
  # 添加汇总工作表 (Add summary worksheet)
  addWorksheet(wb, "Summary")
  
  # 过滤掉Significant为0的数据用于统计
  # Filter out data with Significant=0 for statistics
  filtered_enrichment_data <- enrichment_data %>% filter(Significant > 0)
  
  # 计算每个组每个本体的GO条目数 (Calculate number of GO terms for each group and ontology)
  # 使用pivot_wider替代spread函数
  if(nrow(filtered_enrichment_data) > 0) {
    summary_data <- filtered_enrichment_data %>%
      group_by(Group, ontology) %>%
      summarise(Term_Count = n(), .groups = 'drop') %>%
      pivot_wider(names_from = ontology, values_from = Term_Count, values_fill = 0)
    
    # 如果没有summary数据，添加提示信息 (If no summary data, add information)
    if(nrow(summary_data) == 0) {
      writeData(wb, "Summary", data.frame(Message = "No significant GO terms found"))
    } else {
      # 确保所有本体类型都存在，如果不存在则添加0值列
      # (Ensure all ontology types exist, add 0-value columns if not)
      for(ont in ontologies) {
        if(!(ont %in% colnames(summary_data))) {
          summary_data[[ont]] <- 0
        }
      }
      
      # 写入汇总数据 (Write summary data)
      writeData(wb, "Summary", summary_data)
      setColWidths(wb, "Summary", cols = 1:ncol(summary_data), widths = "auto")
    }
  } else {
    writeData(wb, "Summary", data.frame(Message = "No significant GO terms found after filtering Significant=0"))
  }
  
  # 添加包含所有过滤后结果的工作表 (Add worksheet with all filtered results)
  addWorksheet(wb, "All_Results")
  
  if(nrow(filtered_enrichment_data) > 0) {
    # 移除额外生成的对数列和p值列
    extra_cols <- c("weight_log", "classic_log", "KS_log", "classic_pvalue", "KS_pvalue")
    if(any(extra_cols %in% colnames(filtered_enrichment_data))) {
      filtered_enrichment_data <- filtered_enrichment_data[, !names(filtered_enrichment_data) %in% extra_cols]
    }
    
    # 写入所有数据 (Write all data)
    writeData(wb, "All_Results", filtered_enrichment_data)
    
    # 自动调整列宽 (Auto-adjust column widths)
    setColWidths(wb, "All_Results", cols = 1:ncol(filtered_enrichment_data), widths = "auto")
  } else {
    writeData(wb, "All_Results", data.frame(Message = "No significant GO terms after filtering Significant=0"))
  }
  
  # 保存工作簿 (Save workbook)
  saveWorkbook(wb, output_file, overwrite = TRUE)
  log_message(paste0("结果已保存到: ", output_file, " (Results saved to: ", output_file, ")"))
}

#########################################################################################################################
# 主函数 (Main function)
#########################################################################################################################

main <- function() {
  log_message("开始主函数 (Starting main function)")
  
  # 解析命令行参数 (Parse command line arguments)
  args <- commandArgs(trailingOnly = TRUE)
  
  # 设置输入文件路径 (Set input file paths)
  species_prefix <- "Aty"  # 物种前缀，用于自动构建文件名
  genes_file <- paste0(species_prefix, "_processed_raw_genes_result")
  numbers_file <- paste0(species_prefix, "_processed_raw_numbers_result")
  gene2go_file <- paste0(species_prefix, "_gene2go")
  all_genes_file <- paste0(species_prefix, "_all_genes_list")
  output_dir <- paste0(species_prefix, "_CTDGs_size_enrichment_results")
  
  # 如果提供了命令行参数(If command line arguments are provided)
  if(length(args) >= 1) genes_file <- args[1]
  if(length(args) >= 2) numbers_file <- args[2]
  if(length(args) >= 3) gene2go_file <- args[3]
  if(length(args) >= 4) all_genes_file <- args[4]
  if(length(args) >= 5) output_dir <- args[5]
  
  # 创建输出目录 (Create output directory)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 设置输出文件路径 (Set output file paths)
  output_xlsx <- file.path(output_dir, "CTDGs_size_enrichment_results.xlsx")
  output_mf_plot <- file.path(output_dir, "CTDGs_size_enrichment_MF.pdf")
  output_bp_plot <- file.path(output_dir, "CTDGs_size_enrichment_BP.pdf")
  output_cc_plot <- file.path(output_dir, "CTDGs_size_enrichment_CC.pdf")
  
  # 读取输入文件 (Read input files)
  log_message("读取输入文件 (Reading input files)")
  tryCatch({
    gene2go <- read_gene2go(gene2go_file)
    if(length(gene2go) == 0) {
      log_message("错误: 无法读取gene2go映射文件，终止处理 (Error: Unable to read gene2go mapping file, terminating)")
      return()
    }
    
    genes_data <- read_ctdg_genes(genes_file)
    if(nrow(genes_data) == 0) {
      log_message("错误: 无法读取CTDG基因数据，终止处理 (Error: Unable to read CTDG gene data, terminating)")
      return()
    }
    
    numbers_data <- read_ctdg_numbers(numbers_file)
    if(nrow(numbers_data) == 0) {
      log_message("错误: 无法读取CTDG簇信息，终止处理 (Error: Unable to read CTDG cluster information, terminating)")
      return()
    }
    
    all_genes <- read_all_genes(all_genes_file)
    if(length(all_genes) == 0) {
      log_message("错误: 无法读取所有基因列表，终止处理 (Error: Unable to read all genes list, terminating)")
      return()
    }
    
    # 根据基因簇大小对基因进行分组 (Group genes by cluster size)
    gene_groups <- group_genes_by_cluster_size(genes_data, numbers_data)
    if(length(gene_groups) == 0) {
      log_message("错误: 基因簇大小分组失败，终止处理 (Error: Failed to group genes by cluster size, terminating)")
      return()
    }
    
    # 执行GO富集分析 (Perform GO enrichment analysis)
    log_message("开始执行GO富集分析 (Starting GO enrichment analysis)")
    enrichment_results <- perform_enrichment_for_all_groups(gene_groups, all_genes, gene2go)
    
    # 检查是否有结果 (Check if there are results)
    if(nrow(enrichment_results) > 0) {
      # 创建可视化 (Create visualizations)
      log_message("创建可视化 (Creating visualizations)")
      mf_plot <- create_go_dotplot(enrichment_results, "MF", output_mf_plot)
      bp_plot <- create_go_dotplot(enrichment_results, "BP", output_bp_plot)
      cc_plot <- create_go_dotplot(enrichment_results, "CC", output_cc_plot)
      
      # 保存结果到Excel文件 (Save results to Excel file)
      save_to_excel(enrichment_results, output_xlsx)
      
      log_message(paste0("处理完成，结果已保存到目录: ", output_dir))
    } else {
      log_message("警告: 没有找到任何显著富集的GO条目 (Warning: No significant GO terms found)")
      
      # 创建空的点图 (Create empty plots)
      for(ont in c("MF", "BP", "CC")) {
        output_file <- file.path(output_dir, paste0("CTDGs_size_enrichment_", ont, ".pdf"))
        p <- ggplot() + 
          theme_void() + 
          annotate("text", x = 0, y = 0, label = paste0("No significant ", ont, " terms found")) +
          theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
        ggsave(output_file, plot = p, width = 8.26, height = 9, units = "in")
      }
      
      # 创建空的Excel文件 (Create empty Excel file)
      wb <- createWorkbook()
      addWorksheet(wb, "Info")
      writeData(wb, "Info", data.frame(Message = "No significant GO terms found"))
      saveWorkbook(wb, output_xlsx, overwrite = TRUE)
      
      log_message(paste0("处理完成，但没有找到显著结果，已创建空白输出文件: ", output_dir))
    }
  }, error = function(e) {
    log_message(paste0("错误: 在主函数执行过程中发生错误: ", e$message, " (Error occurred during main function execution)"))
  })
  
  log_message("脚本执行完成 (Script execution completed)")
}

# 执行主函数 (Execute main function)
main()