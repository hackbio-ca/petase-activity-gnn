# Comparison Analysis Script
# This script compares the processed results with the existing data

library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(corrplot)

#' Compare processed results with existing data
#' 
#' @param existing_data_file Path to plate_data_subset.csv
#' @param processed_data_file Path to processed_plate_data.csv
#' @param tetramer_data_file Path to tetramer averages file
#' @return List of comparison results
compare_analysis_results <- function(existing_data_file = "plate_data_subset.csv",
                                    processed_data_file = "processed_plate_data.csv",
                                    tetramer_data_file = "tetramer_averages_1x4.csv") {
  
  cat("=== Comparison Analysis ===\n")
  
  # Read the data files
  cat("Reading data files...\n")
  
  # Existing data
  if (!file.exists(existing_data_file)) {
    stop("Existing data file not found: ", existing_data_file)
  }
  existing_data <- read_csv(existing_data_file, show_col_types = FALSE)
  
  # Processed data
  if (!file.exists(processed_data_file)) {
    stop("Processed data file not found: ", processed_data_file)
  }
  processed_data <- read_csv(processed_data_file, show_col_types = FALSE)
  
  # Tetramer data (optional)
  tetramer_data <- NULL
  if (file.exists(tetramer_data_file)) {
    tetramer_data <- read_csv(tetramer_data_file, show_col_types = FALSE)
  }
  
  cat("✓ Data files loaded successfully\n")
  cat("  - Existing data:", nrow(existing_data), "rows\n")
  cat("  - Processed data:", nrow(processed_data), "rows\n")
  if (!is.null(tetramer_data)) {
    cat("  - Tetramer data:", nrow(tetramer_data), "rows\n")
  }
  
  # Perform comparisons
  comparison_results <- list()
  
  # 1. Data structure comparison
  cat("\n1. Comparing data structures...\n")
  comparison_results$structure <- compare_data_structures(existing_data, processed_data)
  
  # 2. Statistical comparison
  cat("2. Performing statistical comparisons...\n")
  comparison_results$statistics <- compare_statistics(existing_data, processed_data)
  
  # 3. Gene-level comparison
  cat("3. Comparing gene-level results...\n")
  comparison_results$genes <- compare_genes(existing_data, processed_data, tetramer_data)
  
  # 4. Create comparison visualizations
  cat("4. Creating comparison plots...\n")
  create_comparison_plots(existing_data, processed_data, tetramer_data)
  
  # 5. Generate comparison report
  cat("5. Generating comparison report...\n")
  generate_comparison_report(comparison_results, existing_data, processed_data, tetramer_data)
  
  cat("✓ Comparison analysis complete\n")
  
  return(comparison_results)
}

#' Compare data structures between existing and processed data
#' 
#' @param existing_data Existing data
#' @param processed_data Processed data
#' @return List of structure comparison results
compare_data_structures <- function(existing_data, processed_data) {
  
  structure_comparison <- list(
    existing_columns = colnames(existing_data),
    processed_columns = colnames(processed_data),
    common_columns = intersect(colnames(existing_data), colnames(processed_data)),
    existing_only = setdiff(colnames(existing_data), colnames(processed_data)),
    processed_only = setdiff(colnames(processed_data), colnames(existing_data)),
    existing_genes = unique(existing_data$gene),
    processed_genes = unique(processed_data$gene),
    common_genes = intersect(unique(existing_data$gene), unique(processed_data$gene)),
    existing_plates = unique(existing_data$plate),
    processed_plates = unique(processed_data$plate)
  )
  
  cat("  - Common columns:", length(structure_comparison$common_columns), "\n")
  cat("  - Existing only columns:", length(structure_comparison$existing_only), "\n")
  cat("  - Processed only columns:", length(structure_comparison$processed_only), "\n")
  cat("  - Common genes:", length(structure_comparison$common_genes), "\n")
  
  return(structure_comparison)
}

#' Compare statistical distributions
#' 
#' @param existing_data Existing data
#' @param processed_data Processed data
#' @return List of statistical comparison results
compare_statistics <- function(existing_data, processed_data) {
  
  # Compare readout values
  existing_stats <- list(
    mean = mean(existing_data$readout_value, na.rm = TRUE),
    median = median(existing_data$readout_value, na.rm = TRUE),
    sd = sd(existing_data$readout_value, na.rm = TRUE),
    min = min(existing_data$readout_value, na.rm = TRUE),
    max = max(existing_data$readout_value, na.rm = TRUE),
    q25 = quantile(existing_data$readout_value, 0.25, na.rm = TRUE),
    q75 = quantile(existing_data$readout_value, 0.75, na.rm = TRUE)
  )
  
  processed_stats <- list(
    mean = mean(processed_data$readout_value, na.rm = TRUE),
    median = median(processed_data$readout_value, na.rm = TRUE),
    sd = sd(processed_data$readout_value, na.rm = TRUE),
    min = min(processed_data$readout_value, na.rm = TRUE),
    max = max(processed_data$readout_value, na.rm = TRUE),
    q25 = quantile(processed_data$readout_value, 0.25, na.rm = TRUE),
    q75 = quantile(processed_data$readout_value, 0.75, na.rm = TRUE)
  )
  
  # Calculate differences
  stat_differences <- list(
    mean_diff = processed_stats$mean - existing_stats$mean,
    median_diff = processed_stats$median - existing_stats$median,
    sd_ratio = processed_stats$sd / existing_stats$sd,
    range_existing = existing_stats$max - existing_stats$min,
    range_processed = processed_stats$max - processed_stats$min
  )
  
  cat("  - Mean difference:", round(stat_differences$mean_diff, 4), "\n")
  cat("  - Median difference:", round(stat_differences$median_diff, 4), "\n")
  cat("  - SD ratio (processed/existing):", round(stat_differences$sd_ratio, 4), "\n")
  
  return(list(
    existing = existing_stats,
    processed = processed_stats,
    differences = stat_differences
  ))
}

#' Compare gene-level results
#' 
#' @param existing_data Existing data
#' @param processed_data Processed data
#' @param tetramer_data Tetramer data (optional)
#' @return List of gene comparison results
compare_genes <- function(existing_data, processed_data, tetramer_data = NULL) {
  
  # Gene-level summaries for existing data
  existing_gene_stats <- existing_data %>%
    group_by(gene) %>%
    summarise(
      count = n(),
      mean_readout = mean(readout_value, na.rm = TRUE),
      median_readout = median(readout_value, na.rm = TRUE),
      sd_readout = sd(readout_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(data_source = "existing")
  
  # Gene-level summaries for processed data
  processed_gene_stats <- processed_data %>%
    group_by(gene) %>%
    summarise(
      count = n(),
      mean_readout = mean(readout_value, na.rm = TRUE),
      median_readout = median(readout_value, na.rm = TRUE),
      sd_readout = sd(readout_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(data_source = "processed")
  
  # Tetramer gene summaries (if available)
  tetramer_gene_stats <- NULL
  if (!is.null(tetramer_data)) {
    tetramer_gene_stats <- tetramer_data %>%
      group_by(gene) %>%
      summarise(
        count = n(),
        mean_readout = mean(readout_value, na.rm = TRUE),
        median_readout = median(readout_value, na.rm = TRUE),
        sd_readout = sd(readout_value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(data_source = "tetramer")
  }
  
  # Combine for comparison
  combined_gene_stats <- bind_rows(existing_gene_stats, processed_gene_stats, tetramer_gene_stats)
  
  # Find genes present in both datasets
  common_genes <- intersect(existing_gene_stats$gene, processed_gene_stats$gene)
  
  # Calculate correlations for common genes
  gene_correlations <- NULL
  if (length(common_genes) > 1) {
    existing_means <- existing_gene_stats %>% filter(gene %in% common_genes) %>% pull(mean_readout)
    processed_means <- processed_gene_stats %>% filter(gene %in% common_genes) %>% pull(mean_readout)
    
    if (length(existing_means) == length(processed_means)) {
      gene_correlations <- list(
        pearson = cor(existing_means, processed_means, use = "complete.obs"),
        spearman = cor(existing_means, processed_means, method = "spearman", use = "complete.obs")
      )
      cat("  - Gene mean correlation (Pearson):", round(gene_correlations$pearson, 4), "\n")
      cat("  - Gene mean correlation (Spearman):", round(gene_correlations$spearman, 4), "\n")
    }
  }
  
  return(list(
    combined_stats = combined_gene_stats,
    common_genes = common_genes,
    correlations = gene_correlations
  ))
}

#' Create comparison visualization plots
#' 
#' @param existing_data Existing data
#' @param processed_data Processed data
#' @param tetramer_data Tetramer data (optional)
create_comparison_plots <- function(existing_data, processed_data, tetramer_data = NULL) {
  
  # 1. Distribution comparison
  combined_data <- bind_rows(
    existing_data %>% mutate(data_source = "Existing"),
    processed_data %>% mutate(data_source = "Processed")
  )
  
  p1 <- ggplot(combined_data, aes(x = readout_value, fill = data_source)) +
    geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
    facet_wrap(~data_source, ncol = 1) +
    labs(title = "Distribution Comparison: Existing vs Processed",
         x = "Readout Value", y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("comparison_distributions.png", p1, width = 10, height = 8, dpi = 300)
  
  # 2. Gene-level comparison
  gene_comparison_data <- bind_rows(
    existing_data %>% 
      group_by(gene) %>% 
      summarise(mean_readout = mean(readout_value, na.rm = TRUE), .groups = "drop") %>%
      mutate(data_source = "Existing"),
    processed_data %>% 
      group_by(gene) %>% 
      summarise(mean_readout = mean(readout_value, na.rm = TRUE), .groups = "drop") %>%
      mutate(data_source = "Processed")
  )
  
  p2 <- ggplot(gene_comparison_data, aes(x = reorder(gene, mean_readout), 
                                        y = mean_readout, fill = data_source)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(title = "Gene Mean Readout Comparison",
         x = "Gene", y = "Mean Readout Value") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("comparison_genes.png", p2, width = 12, height = 10, dpi = 300)
  
  # 3. Scatter plot for correlation (if data allows)
  common_genes <- intersect(unique(existing_data$gene), unique(processed_data$gene))
  
  if (length(common_genes) > 1) {
    existing_gene_means <- existing_data %>%
      filter(gene %in% common_genes) %>%
      group_by(gene) %>%
      summarise(existing_mean = mean(readout_value, na.rm = TRUE), .groups = "drop")
    
    processed_gene_means <- processed_data %>%
      filter(gene %in% common_genes) %>%
      group_by(gene) %>%
      summarise(processed_mean = mean(readout_value, na.rm = TRUE), .groups = "drop")
    
    correlation_data <- existing_gene_means %>%
      inner_join(processed_gene_means, by = "gene")
    
    p3 <- ggplot(correlation_data, aes(x = existing_mean, y = processed_mean)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
      labs(title = "Gene Mean Readout Correlation",
           x = "Existing Data Mean", y = "Processed Data Mean") +
      theme_minimal()
    
    ggsave("comparison_correlation.png", p3, width = 8, height = 8, dpi = 300)
  }
  
  # 4. Include tetramer comparison if available
  if (!is.null(tetramer_data)) {
    tetramer_gene_means <- tetramer_data %>%
      group_by(gene) %>%
      summarise(tetramer_mean = mean(readout_value, na.rm = TRUE), .groups = "drop")
    
    all_gene_comparison <- existing_data %>%
      group_by(gene) %>%
      summarise(existing_mean = mean(readout_value, na.rm = TRUE), .groups = "drop") %>%
      full_join(processed_data %>%
                  group_by(gene) %>%
                  summarise(processed_mean = mean(readout_value, na.rm = TRUE), .groups = "drop"),
                by = "gene") %>%
      full_join(tetramer_gene_means, by = "gene") %>%
      filter(!is.na(gene))
    
    # Create comparison plot with all three
    comparison_long <- all_gene_comparison %>%
      select(gene, existing_mean, processed_mean, tetramer_mean) %>%
      tidyr::pivot_longer(cols = c(existing_mean, processed_mean, tetramer_mean),
                         names_to = "data_type", values_to = "mean_readout") %>%
      filter(!is.na(mean_readout))
    
    p4 <- ggplot(comparison_long, aes(x = reorder(gene, mean_readout), 
                                     y = mean_readout, fill = data_type)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      labs(title = "Complete Comparison: Existing vs Processed vs Tetramer",
           x = "Gene", y = "Mean Readout Value") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave("comparison_complete.png", p4, width = 14, height = 10, dpi = 300)
  }
  
  cat("✓ Comparison plots saved\n")
}

#' Generate comprehensive comparison report
#' 
#' @param comparison_results Results from comparison analysis
#' @param existing_data Existing data
#' @param processed_data Processed data
#' @param tetramer_data Tetramer data
generate_comparison_report <- function(comparison_results, existing_data, processed_data, tetramer_data = NULL) {
  
  report_content <- c(
    "# Halo Assay Analysis Comparison Report",
    paste("Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Data Overview",
    paste("- Existing data rows:", nrow(existing_data)),
    paste("- Processed data rows:", nrow(processed_data)),
    if (!is.null(tetramer_data)) paste("- Tetramer data rows:", nrow(tetramer_data)) else NULL,
    "",
    "## Structure Comparison",
    paste("- Common columns:", length(comparison_results$structure$common_columns)),
    paste("- Columns only in existing:", paste(comparison_results$structure$existing_only, collapse = ", ")),
    paste("- Columns only in processed:", paste(comparison_results$structure$processed_only, collapse = ", ")),
    paste("- Common genes:", length(comparison_results$structure$common_genes)),
    "",
    "## Statistical Comparison",
    "### Readout Values",
    paste("- Mean difference (processed - existing):", 
          round(comparison_results$statistics$differences$mean_diff, 6)),
    paste("- Median difference (processed - existing):", 
          round(comparison_results$statistics$differences$median_diff, 6)),
    paste("- Standard deviation ratio (processed/existing):", 
          round(comparison_results$statistics$differences$sd_ratio, 4)),
    "",
    "### Existing Data Statistics",
    paste("- Mean:", round(comparison_results$statistics$existing$mean, 6)),
    paste("- Median:", round(comparison_results$statistics$existing$median, 6)),
    paste("- Standard deviation:", round(comparison_results$statistics$existing$sd, 6)),
    paste("- Range:", round(comparison_results$statistics$existing$min, 6), "to", 
          round(comparison_results$statistics$existing$max, 6)),
    "",
    "### Processed Data Statistics",
    paste("- Mean:", round(comparison_results$statistics$processed$mean, 6)),
    paste("- Median:", round(comparison_results$statistics$processed$median, 6)),
    paste("- Standard deviation:", round(comparison_results$statistics$processed$sd, 6)),
    paste("- Range:", round(comparison_results$statistics$processed$min, 6), "to", 
          round(comparison_results$statistics$processed$max, 6)),
    ""
  )
  
  # Add gene correlation information if available
  if (!is.null(comparison_results$genes$correlations)) {
    report_content <- c(report_content,
      "## Gene-Level Correlation",
      paste("- Pearson correlation:", round(comparison_results$genes$correlations$pearson, 4)),
      paste("- Spearman correlation:", round(comparison_results$genes$correlations$spearman, 4)),
      ""
    )
  }
  
  # Add top performing genes comparison
  top_existing <- comparison_results$genes$combined_stats %>%
    filter(data_source == "existing") %>%
    arrange(desc(mean_readout)) %>%
    head(5)
  
  top_processed <- comparison_results$genes$combined_stats %>%
    filter(data_source == "processed") %>%
    arrange(desc(mean_readout)) %>%
    head(5)
  
  report_content <- c(report_content,
    "## Top Performing Genes",
    "### Existing Data (Top 5)",
    paste("1.", top_existing$gene[1], "- Mean:", round(top_existing$mean_readout[1], 4)),
    paste("2.", top_existing$gene[2], "- Mean:", round(top_existing$mean_readout[2], 4)),
    paste("3.", top_existing$gene[3], "- Mean:", round(top_existing$mean_readout[3], 4)),
    paste("4.", top_existing$gene[4], "- Mean:", round(top_existing$mean_readout[4], 4)),
    paste("5.", top_existing$gene[5], "- Mean:", round(top_existing$mean_readout[5], 4)),
    "",
    "### Processed Data (Top 5)",
    paste("1.", top_processed$gene[1], "- Mean:", round(top_processed$mean_readout[1], 4)),
    paste("2.", top_processed$gene[2], "- Mean:", round(top_processed$mean_readout[2], 4)),
    paste("3.", top_processed$gene[3], "- Mean:", round(top_processed$mean_readout[3], 4)),
    paste("4.", top_processed$gene[4], "- Mean:", round(top_processed$mean_readout[4], 4)),
    paste("5.", top_processed$gene[5], "- Mean:", round(top_processed$mean_readout[5], 4)),
    "",
    "## Files Generated",
    "- comparison_distributions.png: Distribution comparison plot",
    "- comparison_genes.png: Gene-level comparison plot",
    "- comparison_correlation.png: Gene correlation scatter plot",
    if (!is.null(tetramer_data)) "- comparison_complete.png: Complete three-way comparison" else NULL,
    "",
    "## Methodology",
    "1. Structural comparison of data formats and content",
    "2. Statistical comparison of readout value distributions",
    "3. Gene-level performance comparison and correlation analysis",
    "4. Visual comparison through multiple plot types",
    "5. Tetramer averaging validation (if applicable)"
  )
  
  # Write the report
  writeLines(report_content, "comparison_report.md")
  cat("✓ Comparison report saved to comparison_report.md\n")
}

# If running this script directly
if (!interactive()) {
  cat("Running comparison analysis...\n")
  results <- compare_analysis_results()
  cat("Comparison analysis complete!\n")
}