# Batch Processing Script for Halo Assay Analysis
# This script orchestrates the complete analysis pipeline

library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)

# Source the required scripts
source("plate_analysis_main.R")
source("tetramer_averaging.R")

#' Complete batch processing pipeline
#' 
#' This function runs the complete analysis pipeline:
#' 1. Process all plate images
#' 2. Calculate tetramer averages 
#' 3. Generate summary reports
#' 
#' @param base_dir Base directory containing the analysis files
#' @param tetramer_method Method for tetramer calculation ("2x2" - standard layout)
#' @param results_dir Directory to save results (default: current directory)
#' @return List containing all analysis results
run_complete_analysis <- function(base_dir = ".", tetramer_method = "2x2", results_dir = ".") {
  
  cat("=== Starting Complete Halo Assay Analysis ===\n")
  cat("Base directory:", base_dir, "\n")
  cat("Tetramer method:", tetramer_method, "\n")
  cat("Results directory:", results_dir, "\n")
  cat("Analysis started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Step 1: Process all plate images
  cat("STEP 1: Processing all plate images...\n")
  tryCatch({
    raw_results <- process_all_plates(base_dir)
    cat("✓ Successfully processed", nrow(raw_results), "colonies\n\n")
  }, error = function(e) {
    cat("✗ Error in plate processing:", e$message, "\n")
    stop("Cannot continue without plate processing results")
  })
  
  # Step 2: Calculate tetramer averages
  cat("STEP 2: Calculating tetramer averages...\n")
  tryCatch({
    tetramer_results <- process_tetramer_averages(
      input_file = "processed_plate_data.csv",
      method = tetramer_method
    )
    cat("✓ Successfully calculated", nrow(tetramer_results), "tetramer averages\n\n")
  }, error = function(e) {
    cat("✗ Error in tetramer calculation:", e$message, "\n")
    tetramer_results <- NULL
  })
  
  # Step 3: Generate summary statistics
  cat("STEP 3: Generating summary statistics...\n")
  summary_stats <- generate_summary_statistics(raw_results, tetramer_results)
  
  # Step 4: Create visualizations
  cat("STEP 4: Creating visualizations...\n")
  create_analysis_plots(raw_results, tetramer_results)
  
  # Step 5: Generate final report
  cat("STEP 5: Generating analysis report...\n")
  generate_analysis_report(raw_results, tetramer_results, summary_stats)
  
  cat("=== Analysis Complete ===\n")
  cat("Analysis finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  return(list(
    raw_data = raw_results,
    tetramer_data = tetramer_results,
    summary = summary_stats
  ))
}

#' Generate comprehensive summary statistics
#' 
#' @param raw_data Raw colony data
#' @param tetramer_data Tetramer average data
#' @return List of summary statistics
generate_summary_statistics <- function(raw_data, tetramer_data = NULL) {
  
  cat("Generating summary statistics...\n")
  
  # Overall statistics for raw data
  raw_summary <- list(
    total_colonies = nrow(raw_data),
    unique_plates = length(unique(raw_data$plate)),
    unique_genes = length(unique(raw_data$gene)),
    readout_range = range(raw_data$readout_value, na.rm = TRUE),
    readout_mean = mean(raw_data$readout_value, na.rm = TRUE),
    readout_median = median(raw_data$readout_value, na.rm = TRUE),
    colony_size_range = range(raw_data$colony_size, na.rm = TRUE),
    colony_size_mean = mean(raw_data$colony_size, na.rm = TRUE)
  )
  
  # Gene-level statistics
  gene_stats <- raw_data %>%
    group_by(gene) %>%
    summarise(
      count = n(),
      mean_readout = mean(readout_value, na.rm = TRUE),
      median_readout = median(readout_value, na.rm = TRUE),
      sd_readout = sd(readout_value, na.rm = TRUE),
      mean_colony_size = mean(colony_size, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_readout))
  
  # Plate-level statistics
  plate_stats <- raw_data %>%
    group_by(plate) %>%
    summarise(
      count = n(),
      mean_readout = mean(readout_value, na.rm = TRUE),
      median_readout = median(readout_value, na.rm = TRUE),
      sd_readout = sd(readout_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Tetramer statistics (if available)
  tetramer_summary <- NULL
  if (!is.null(tetramer_data)) {
    tetramer_summary <- list(
      total_tetramers = nrow(tetramer_data),
      readout_range = range(tetramer_data$readout_value, na.rm = TRUE),
      readout_mean = mean(tetramer_data$readout_value, na.rm = TRUE),
      readout_median = median(tetramer_data$readout_value, na.rm = TRUE)
    )
  }
  
  # Save summary statistics
  summary_data <- list(
    raw_summary = raw_summary,
    gene_stats = gene_stats,
    plate_stats = plate_stats,
    tetramer_summary = tetramer_summary
  )
  
  # Write gene and plate statistics to files
  write_csv(gene_stats, "summary_by_gene.csv")
  write_csv(plate_stats, "summary_by_plate.csv")
  
  cat("✓ Summary statistics generated and saved\n")
  
  return(summary_data)
}

#' Create analysis plots and visualizations
#' 
#' @param raw_data Raw colony data
#' @param tetramer_data Tetramer average data
create_analysis_plots <- function(raw_data, tetramer_data = NULL) {
  
  cat("Creating visualization plots...\n")
  
  # 1. Distribution of readout values
  p1 <- ggplot(raw_data, aes(x = readout_value)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    labs(title = "Distribution of Readout Values (Raw Data)",
         x = "Readout Value", y = "Count") +
    theme_minimal()
  
  ggsave("plot_readout_distribution.png", p1, width = 10, height = 6, dpi = 300)
  
  # 2. Readout values by gene
  p2 <- ggplot(raw_data, aes(x = reorder(gene, readout_value, FUN = median), 
                             y = readout_value)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Readout Values by Gene",
         x = "Gene", y = "Readout Value") +
    theme_minimal()
  
  ggsave("plot_readout_by_gene.png", p2, width = 12, height = 8, dpi = 300)
  
  # 3. Colony size vs readout value
  p3 <- ggplot(raw_data, aes(x = colony_size, y = readout_value)) +
    geom_point(alpha = 0.5, color = "darkblue") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(title = "Colony Size vs Readout Value",
         x = "Colony Size", y = "Readout Value") +
    theme_minimal()
  
  ggsave("plot_size_vs_readout.png", p3, width = 10, height = 6, dpi = 300)
  
  # 4. Plate-to-plate comparison
  p4 <- ggplot(raw_data, aes(x = plate, y = readout_value)) +
    geom_boxplot(fill = "orange", alpha = 0.7) +
    labs(title = "Readout Values by Plate",
         x = "Plate", y = "Readout Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("plot_readout_by_plate.png", p4, width = 10, height = 6, dpi = 300)
  
  # 5. Tetramer comparison (if available)
  if (!is.null(tetramer_data)) {
    p5 <- ggplot(tetramer_data, aes(x = reorder(gene, readout_value, FUN = median), 
                                   y = readout_value)) +
      geom_boxplot(fill = "lightgreen", alpha = 0.7) +
      coord_flip() +
      labs(title = "Tetramer Average Readout Values by Gene",
           x = "Gene", y = "Tetramer Average Readout Value") +
      theme_minimal()
    
    ggsave("plot_tetramer_by_gene.png", p5, width = 12, height = 8, dpi = 300)
  }
  
  cat("✓ Visualization plots created and saved\n")
}

#' Generate comprehensive analysis report
#' 
#' @param raw_data Raw colony data
#' @param tetramer_data Tetramer average data
#' @param summary_stats Summary statistics
generate_analysis_report <- function(raw_data, tetramer_data, summary_stats) {
  
  cat("Generating analysis report...\n")
  
  report_content <- c(
    "# Halo Assay Analysis Report",
    paste("Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Summary",
    paste("- Total colonies analyzed:", summary_stats$raw_summary$total_colonies),
    paste("- Number of plates:", summary_stats$raw_summary$unique_plates),
    paste("- Number of unique genes:", summary_stats$raw_summary$unique_genes),
    paste("- Readout value range:", paste(round(summary_stats$raw_summary$readout_range, 4), collapse = " to ")),
    paste("- Mean readout value:", round(summary_stats$raw_summary$readout_mean, 4)),
    paste("- Median readout value:", round(summary_stats$raw_summary$readout_median, 4)),
    "",
    "## Top Performing Genes (by mean readout)",
    ""
  )
  
  # Add top genes table
  top_genes <- head(summary_stats$gene_stats, 10)
  for (i in 1:nrow(top_genes)) {
    report_content <- c(report_content,
      paste(i, ".", top_genes$gene[i], "- Mean:", round(top_genes$mean_readout[i], 4),
            "Count:", top_genes$count[i]))
  }
  
  # Add tetramer results if available
  if (!is.null(tetramer_data)) {
    report_content <- c(report_content,
      "",
      "## Tetramer Analysis Results",
      paste("- Total tetramers:", summary_stats$tetramer_summary$total_tetramers),
      paste("- Tetramer readout range:", paste(round(summary_stats$tetramer_summary$readout_range, 4), collapse = " to ")),
      paste("- Mean tetramer readout:", round(summary_stats$tetramer_summary$readout_mean, 4)),
      paste("- Median tetramer readout:", round(summary_stats$tetramer_summary$readout_median, 4))
    )
  }
  
  report_content <- c(report_content,
    "",
    "## Files Generated",
    "- processed_plate_data.csv: Raw colony analysis data",
    "- tetramer_averages_1x4.csv: Tetramer averages (if calculated)",
    "- summary_by_gene.csv: Gene-level summary statistics",
    "- summary_by_plate.csv: Plate-level summary statistics",
    "- plot_*.png: Various visualization plots",
    "",
    "## Analysis Method",
    "1. Images processed using HAIP package functions",
    "2. Background correction applied when background images available",
    "3. Tetramer averages calculated from bottom two values in each group",
    "4. Statistical summaries and visualizations generated"
  )
  
  # Write report
  writeLines(report_content, "analysis_report.md")
  cat("✓ Analysis report saved to analysis_report.md\n")
}

# If running this script directly, execute the complete analysis
if (!interactive()) {
  cat("Running complete batch analysis...\n")
  results <- run_complete_analysis()
  cat("Batch analysis complete!\n")
}