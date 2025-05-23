#!/usr/bin/env Rscript

#############################################################
# CTDG-FoldChange Visualization Script (Optimized Version)
# Author: Claude
# Date: 2025-05-17
#
# Creates visualizations for CTDG-FoldChange results.
# This version focuses only on combined plots across conditions.
#
# Main features:
# 1. Combined density plots: Compare fold change distributions between DE-CTDGs and random CTDGs
# 2. Combined boxplots: Visual comparison of expression distributions
#
# Usage:
#   Rscript visualize_ctdg_foldchange.R \
#     --data_dir /path/to/analysis_output \
#     --output_dir /path/to/output \
#     --conditions "condition1,condition2" \
#     --size_groups "member=2,3<member<5,6<member<10" \
#     --color_palette "viridis"
#############################################################

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(jsonlite)
  library(optparse)
  library(viridis)  # For color palettes
})

# Parse command line arguments
parse_arguments <- function() {
  option_list <- list(
    make_option("--data_dir", type="character", default=NULL, 
                help="Data directory containing Python analysis script output"),
    make_option("--output_dir", type="character", default="./plots", 
                help="Output directory for generated plots"),
    make_option("--config_file", type="character", default=NULL,
                help="Config file path (JSON format), overrides other parameters if provided"),
    make_option("--conditions", type="character", default=NULL,
                help="Conditions to process, comma-separated"),
    make_option("--size_groups", type="character", default="member=2,3<member<5,6<member<10,member>10",
                help="Size groups to process, comma-separated"),
    make_option("--color_palette", type="character", default="viridis",
                help="Color palette, options: viridis, brewer, grey"),
    make_option("--width", type="integer", default=8.3,
                help="Plot width (inches)"),
    make_option("--height", type="integer", default=4,
                help="Plot height (inches)")
  )
  
  parser <- OptionParser(option_list=option_list,
                         description="Visualize CTDG-FoldChange analysis results")
  
  args <- parse_args(parser)
  
  # Load parameters from config file if provided
  if (!is.null(args$config_file)) {
    if (file.exists(args$config_file)) {
      config <- fromJSON(args$config_file)
      
      # Update parameters
      if ("data_dir" %in% names(config)) args$data_dir <- config$data_dir
      if ("conditions" %in% names(config)) args$conditions <- paste(config$conditions, collapse=",")
      if ("size_groups" %in% names(config)) args$size_groups <- paste(config$size_groups, collapse=",")
    } else {
      stop("Specified config file does not exist: ", args$config_file)
    }
  }
  
  # Check required parameters
  if (is.null(args$data_dir)) {
    stop("Data directory must be provided (--data_dir)")
  }
  
  # Infer conditions from directory if not specified
  if (is.null(args$conditions)) {
    viz_dir <- file.path(args$data_dir, "visualization_data")
    if (dir.exists(viz_dir)) {
      conditions <- list.dirs(viz_dir, full.names=FALSE, recursive=FALSE)
      if (length(conditions) > 0) {
        args$conditions <- paste(conditions, collapse=",")
        cat("Inferred conditions from directory:", args$conditions, "\n")
      } else {
        stop("Unable to infer conditions from directory, please specify with --conditions")
      }
    } else {
      stop("Unable to infer conditions from directory, please specify with --conditions")
    }
  }
  
  # Create output directory
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive=TRUE)
  }
  
  # Parse conditions and size groups
  args$conditions <- strsplit(args$conditions, ",")[[1]]
  args$size_groups <- strsplit(args$size_groups, ",")[[1]]
  
  return(args)
}

# Load visualization data
load_visualization_data <- function(data_dir, condition, verbose=TRUE) {
  if (verbose) cat("Loading data for condition", condition, "\n")
  
  viz_dir <- file.path(data_dir, "visualization_data", condition)
  
  if (!dir.exists(viz_dir)) {
    warning("Directory does not exist: ", viz_dir)
    return(NULL)
  }
  
  # Load density plot data
  density_file <- file.path(viz_dir, "density_data.csv")
  if (file.exists(density_file)) {
    density_data <- read_csv(density_file, show_col_types=FALSE)
    if (verbose) cat("  Loaded density data:", nrow(density_data), "rows\n")
  } else {
    warning("Density data file not found: ", density_file)
    density_data <- NULL
  }
  
  # Load boxplot data
  boxplot_file <- file.path(viz_dir, "boxplot_data.csv")
  if (file.exists(boxplot_file)) {
    boxplot_data <- read_csv(boxplot_file, show_col_types=FALSE)
    if (verbose) cat("  Loaded boxplot data:", nrow(boxplot_data), "rows\n")
  } else {
    warning("Boxplot data file not found: ", boxplot_file)
    boxplot_data <- NULL
  }
  
  # Load summary data
  summary_file <- file.path(viz_dir, "summary_data.json")
  if (file.exists(summary_file)) {
    summary_data <- fromJSON(summary_file)
    if (verbose) cat("  Loaded summary data:", length(summary_data), "size groups\n")
  } else {
    warning("Summary data file not found: ", summary_file)
    summary_data <- NULL
  }
  
  return(list(
    density_data = density_data,
    boxplot_data = boxplot_data,
    summary_data = summary_data
  ))
}

# Create summary table (needed for p-value annotations in density plot)
create_summary_table <- function(summary_data, size_groups, condition) {
  # Ensure data contains specified size groups
  if (is.null(summary_data)) {
    warning("Cannot create summary table for condition", condition, ": Empty data")
    return(NULL)
  }
  
  # Filter specified size groups
  size_groups <- intersect(size_groups, names(summary_data))
  
  if (length(size_groups) == 0) {
    warning("Cannot create summary table for condition", condition, ": No matching size groups")
    return(NULL)
  }
  
  # Create summary table
  summary_rows <- list()
  
  for (size_group in size_groups) {
    if (!size_group %in% names(summary_data)) {
      next
    }
    
    group_data <- summary_data[[size_group]]
    
    # Check if required fields exist
    required_fields <- c(
      "de_ctdg_count", "random_ctdg_count", 
      "de_statistics", "random_statistics", 
      "ks_test", "effect_size"
    )
    
    if (!all(required_fields %in% names(group_data))) {
      missing_fields <- required_fields[!required_fields %in% names(group_data)]
      warning("Size group ", size_group, " missing required fields: ", 
              paste(missing_fields, collapse=", "))
      next
    }
    
    # Check if statistics fields are complete
    if (!all(c("mean", "median") %in% names(group_data$de_statistics)) || 
        !all(c("mean", "median") %in% names(group_data$random_statistics))) {
      warning("Size group ", size_group, " has incomplete statistics data")
      next
    }
    
    # Check if test results are complete
    if (!all(c("statistic", "pvalue") %in% names(group_data$ks_test))) {
      warning("Size group ", size_group, " has incomplete KS test results")
      next
    }
    
    # Check if effect size is complete
    if (!"cohens_d" %in% names(group_data$effect_size)) {
      warning("Size group ", size_group, " missing effect size data")
      next
    }
    
    # Better labels for size groups
    if (size_group == "member=2") {
      size_label <- "Members = 2"
    } else if (size_group == "3<member<5") {
      size_label <- "Members 3-5"
    } else if (size_group == "6<member<10") {
      size_label <- "Members 6-10"
    } else {
      size_label <- "Members >10"
    }
    
    # Extract KS test results
    ks_statistic <- group_data$ks_test$statistic
    ks_pvalue <- group_data$ks_test$pvalue
    
    # Extract effect size
    cohens_d <- group_data$effect_size$cohens_d
    
    # Ensure all values are valid
    if (is.null(group_data$de_ctdg_count) || 
        is.null(group_data$random_ctdg_count) || 
        is.null(group_data$de_statistics$mean) || 
        is.null(group_data$random_statistics$mean) ||
        is.null(group_data$de_statistics$median) || 
        is.null(group_data$random_statistics$median) ||
        is.null(ks_statistic) || 
        is.null(ks_pvalue) || 
        is.null(cohens_d)) {
      warning("Size group ", size_group, " contains NULL values, skipping")
      next
    }
    
    # Create summary row
    tryCatch({
      row <- data.frame(
        size_group = size_label,
        ks_pvalue = ks_pvalue,
        significance = case_when(
          ks_pvalue < 0.001 ~ "***",
          ks_pvalue < 0.01 ~ "**",
          ks_pvalue < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      )
      summary_rows[[size_group]] <- row
    }, error = function(e) {
      warning("Error creating summary row for size group ", size_group, ": ", e$message)
    })
  }
  
  # Merge all rows
  if (length(summary_rows) > 0) {
    # Check if all rows are valid before merging
    valid_rows <- TRUE
    for (size_group in names(summary_rows)) {
      row <- summary_rows[[size_group]]
      if (is.null(row) || nrow(row) == 0) {
        warning("Invalid summary row for size group ", size_group, ", skipping")
        summary_rows[[size_group]] <- NULL
        valid_rows <- FALSE
      }
    }
    
    if (length(summary_rows) == 0) {
      warning("No valid summary rows, returning NULL")
      return(NULL)
    }
    
    # Try to merge all rows
    tryCatch({
      summary_table <- do.call(rbind, summary_rows)
      
      # Sort by size group
      size_order <- c("Members = 2", "Members 3-5", "Members 6-10", "Members >10")
      summary_table$size_group <- factor(summary_table$size_group, levels=size_order)
      summary_table <- summary_table[order(summary_table$size_group),]
      
      return(summary_table)
    }, error = function(e) {
      warning("Error merging summary rows: ", e$message)
      return(NULL)
    })
  } else {
    warning("No summary rows generated, returning NULL")
    return(NULL)
  }
}

# Combine all conditions into one density plot
create_combined_density_plot <- function(all_density_data, all_summary_tables, size_groups, conditions, color_palette) {
  # Check if we have data
  if (length(all_density_data) == 0) {
    warning("No density data available to create combined plot")
    return(NULL)
  }
  
  # Combine all data into one dataframe
  combined_data <- data.frame()
  
  for (condition in names(all_density_data)) {
    data <- all_density_data[[condition]]
    if (!is.null(data) && nrow(data) > 0) {
      # Filter data for specified size groups
      data <- data %>% filter(size_group %in% size_groups)
      
      # Add condition column
      data$condition <- condition
      
      # Append to combined data
      combined_data <- rbind(combined_data, data)
    }
  }
  
  if (nrow(combined_data) == 0) {
    warning("No valid data after filtering and combining")
    return(NULL)
  }
  
  # Set proper order for size groups
  size_group_order <- c("member=2", "3<member<5", "6<member<10", "member>10")
  combined_data$size_group <- factor(combined_data$size_group, levels=size_group_order)
  
  # Set proper order for conditions
  combined_data$condition <- factor(combined_data$condition, levels=conditions)
  
  # Create labels for size groups
  size_group_labels <- c(
    "member=2" = "Members = 2",
    "3<member<5" = "Members 3-5",
    "6<member<10" = "Members 6-10",
    "member>10" = "Members >10"
  )
  
  # Create combined plot
  p <- ggplot(combined_data, aes(x=log2_fold_change, fill=group)) +
    geom_density(alpha=0.7,size=.1) +
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
      axis.text = element_text(size=8,color = "black"),
      axis.line.x = element_line(linewidth = .3,color = "black"),
      axis.line.y = element_line(linewidth = .3,color = "black"),
      legend.title = element_text(size=12),
      legend.text = element_text(size=10),
      strip.text = element_text(size=11, face="bold"),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    geom_vline(xintercept=0, linetype="dashed", color="black",linewidth=.3)
  
  # Add p-value annotations if summary tables are available
  if (!is.null(all_summary_tables) && length(all_summary_tables) > 0) {
    annotations <- data.frame()
    
    for (condition in names(all_summary_tables)) {
      summary_table <- all_summary_tables[[condition]]
      if (!is.null(summary_table) && nrow(summary_table) > 0) {
        # Extract p-values and significance
        condition_anno <- summary_table %>%
          select(size_group, ks_pvalue, significance) %>%
          mutate(
            condition = condition,
            y_pos = 0.95,
            label = paste0("p = ", sprintf("%.3f", ks_pvalue), " ", significance)
          )
        
        # Map size group labels
        size_group_labels_rev <- c(
          "Members = 2" = "member=2",
          "Members 3-5" = "3<member<5",
          "Members 6-10" = "6<member<10",
          "Members >10" = "member>10"
        )
        
        condition_anno$size_group <- size_group_labels_rev[condition_anno$size_group]
        
        # Append to annotations dataframe
        annotations <- rbind(annotations, condition_anno)
      }
    }
    
    if (nrow(annotations) > 0) {
      # Filter out missing values
      annotations <- annotations %>%
        filter(!is.na(size_group)) %>%
        filter(!is.na(condition))
      
      # Set factors for proper ordering
      annotations$size_group <- factor(annotations$size_group, levels=size_group_order)
      annotations$condition <- factor(annotations$condition, levels=conditions)
      
      # Add annotations to plot
      tryCatch({
        p <- p + 
          geom_text(
            data = annotations,
            aes(x = 0, y = y_pos, label = label),
            inherit.aes = FALSE,
            hjust = 0.5,
            size = 3
          )
      }, error = function(e) {
        warning("Error adding annotations to combined plot: ", e$message)
      })
    }
  }
  
  return(p)
}

# Combine all conditions into one boxplot
create_combined_boxplot <- function(all_boxplot_data, size_groups, conditions, color_palette) {
  # Check if we have data
  if (length(all_boxplot_data) == 0) {
    warning("No boxplot data available to create combined plot")
    return(NULL)
  }
  
  # Combine all data into one dataframe
  combined_data <- data.frame()
  
  for (condition in names(all_boxplot_data)) {
    data <- all_boxplot_data[[condition]]
    if (!is.null(data) && nrow(data) > 0) {
      # Filter data for specified size groups
      data <- data %>% filter(size_group %in% size_groups)
      
      # Add condition column
      data$condition <- condition
      
      # Append to combined data
      combined_data <- rbind(combined_data, data)
    }
  }
  
  if (nrow(combined_data) == 0) {
    warning("No valid data after filtering and combining")
    return(NULL)
  }
  
  # Set proper order for size groups
  size_group_order <- c("member=2", "3<member<5", "6<member<10", "member>10")
  combined_data$size_group <- factor(combined_data$size_group, levels=size_group_order)
  
  # Set proper order for conditions
  combined_data$condition <- factor(combined_data$condition, levels=conditions)
  
  # Create labels for size groups
  size_group_labels <- c(
    "member=2" = "Members = 2",
    "3<member<5" = "Members 3-5",
    "6<member<10" = "Members 6-10",
    "member>10" = "Members >10"
  )
  
  # Create combined boxplot
  p <- ggplot(combined_data, aes(x=group, y=log2_fold_change, color=group)) +
    geom_boxplot(alpha=0.8, outlier.size=.3, width=0.6,fill=NA,linewidth=.3) +
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
      axis.text.x = element_text(size=10,color = "black"),
      axis.text.y = element_text(size=7,color = "black"),
      axis.line.x = element_line(linewidth = .3,color = "black"),
      axis.line.y = element_line(linewidth = .3,color = "black"),
      legend.title = element_text(size=12),
      legend.text = element_text(size=10),
      strip.text = element_text(size=11, face="bold"),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    geom_hline(yintercept=0, linetype="dashed", color="black",linewidth=.3)
  
  return(p)
}

# Save plot (PDF only)
save_plot <- function(plot, filename, width, height) {
  # Save PDF version only
  tryCatch({
    pdf(filename, width=width, height=height)
    print(plot)
    dev.off()
    
    cat("Plot saved to:", filename, "\n")
  }, error = function(e) {
    warning("Error saving plot: ", e$message)
  })
}

# Main function
main <- function() {
  # Parse command line arguments
  args <- parse_arguments()
  
  # Try to load config from r_config.json
  config_file <- file.path(args$data_dir, "r_config.json")
  if (file.exists(config_file)) {
    cat("Loading parameters from config file:", config_file, "\n")
    args$config_file <- config_file
    args <- parse_arguments()
  }
  
  cat("Data directory:", args$data_dir, "\n")
  cat("Output directory:", args$output_dir, "\n")
  cat("Processing conditions:", paste(args$conditions, collapse=", "), "\n")
  cat("Processing size groups:", paste(args$size_groups, collapse=", "), "\n")
  
  # Create output directory
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each condition just to collect data
  all_summaries <- list()
  all_density_data <- list()
  all_boxplot_data <- list()
  
  for (condition in args$conditions) {
    cat("\nProcessing condition:", condition, "\n")
    
    # Load data
    viz_data <- load_visualization_data(args$data_dir, condition)
    
    if (is.null(viz_data)) {
      warning("Skipping condition ", condition, ": Unable to load data")
      next
    }
    
    # Store data for combined plots
    if (!is.null(viz_data$density_data) && nrow(viz_data$density_data) > 0) {
      all_density_data[[condition]] <- viz_data$density_data
    }
    
    if (!is.null(viz_data$boxplot_data) && nrow(viz_data$boxplot_data) > 0) {
      all_boxplot_data[[condition]] <- viz_data$boxplot_data
    }
    
    # Create summary table for p-value annotations
    tryCatch({
      summary_table <- create_summary_table(viz_data$summary_data, args$size_groups, condition)
      
      if (!is.null(summary_table) && nrow(summary_table) > 0) {
        # Store for combined plots
        all_summaries[[condition]] <- summary_table
      } else {
        warning("Condition ", condition, " has no valid summary table")
      }
    }, error = function(e) {
      warning("Error creating summary table for condition ", condition, ": ", e$message)
    })
  }
  
  # Create combined density plot
  if (length(all_density_data) > 0) {
    tryCatch({
      combined_density <- create_combined_density_plot(
        all_density_data,
        all_summaries,
        args$size_groups,
        args$conditions,
        args$color_palette
      )
      
      if (!is.null(combined_density)) {
        # Save combined density plot
        combined_density_file <- file.path(args$output_dir, "combined_density_plot.pdf")
        save_plot(combined_density, combined_density_file, args$width, args$height)
      }
    }, error = function(e) {
      warning("Error creating combined density plot: ", e$message)
    })
  }
  
  # Create combined boxplot
  if (length(all_boxplot_data) > 0) {
    tryCatch({
      combined_boxplot <- create_combined_boxplot(
        all_boxplot_data,
        args$size_groups,
        args$conditions,
        args$color_palette
      )
      
      if (!is.null(combined_boxplot)) {
        # Save combined boxplot
        combined_boxplot_file <- file.path(args$output_dir, "combined_boxplot_plot.pdf")
        save_plot(combined_boxplot, combined_boxplot_file, args$width, args$height)
      }
    }, error = function(e) {
      warning("Error creating combined boxplot: ", e$message)
    })
  }
  
  cat("\nAll conditions processed\n")
}

# Execute main function
main()