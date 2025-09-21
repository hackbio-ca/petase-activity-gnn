# Tetramer Averaging Script - 2x2 Grid Layout
# This script calculates the average of the bottom two values for each 2x2 tetramer group
# The plates are laid out in 2x2 tetramers, which is the standard configuration.

library(dplyr)
library(readr)

#' Calculate average of bottom two values for each 2x2 tetramer
#' 
#' For 2x2 tetramer grids, this function:
#' 1. Groups wells into 2x2 grids based on their row and column positions
#' 2. For each 2x2 grid, finds the bottom two values (lowest readout_values)
#' 3. Calculates their average
#' 4. Returns this average as the representative value for the tetramer
#' 
#' @param data Data frame containing plate analysis results
#' @return Data frame with tetramer averages
calculate_tetramer_averages_2x2 <- function(data) {
  
  cat("Calculating tetramer averages for 2x2 grid pattern...\n")
  
  tetramer_results <- data %>%
    # Create 2x2 tetramer groups
    mutate(
      tetramer_col_group = ceiling(column / 2),
      tetramer_row_group = ceiling(row / 2)
    ) %>%
    group_by(plate, plasmid, tetramer_row_group, tetramer_col_group, gene) %>%
    
    # For each 2x2 tetramer, get the bottom two values
    arrange(readout_value) %>%
    slice_head(n = 2) %>%  # Get the 2 lowest values
    
    # Calculate average
    summarise(
      tetramer_avg_readout = mean(readout_value, na.rm = TRUE),
      tetramer_avg_colony_size = mean(colony_size, na.rm = TRUE),
      num_colonies_averaged = n(),
      min_readout = min(readout_value, na.rm = TRUE),
      max_readout = max(readout_value, na.rm = TRUE),
      date_entered = first(date_entered),
      .groups = "drop"
    ) %>%
    
    # Calculate representative positions
    mutate(
      representative_column = (tetramer_col_group - 1) * 2 + 1.5,
      representative_row = (tetramer_row_group - 1) * 2 + 1.5,
      normalization_method = "tetramer_2x2_avg_bottom2"
    ) %>%
    
    # Generate new IDs
    mutate(id = row_number() + 60000) %>%
    
    # Select columns
    select(
      id,
      gene,
      plate,
      plasmid,
      column = representative_column,
      row = representative_row,
      normalization_method,
      readout_value = tetramer_avg_readout,
      colony_size = tetramer_avg_colony_size,
      date_entered,
      tetramer_col_group,
      tetramer_row_group,
      num_colonies_averaged,
      min_readout,
      max_readout
    )
  
  return(tetramer_results)
}

#' Process tetramer averages from input data
#' Process tetramer averages using 2x2 grid method
#' 
#' This function processes the plate data to calculate tetramer averages.
#' The plates are laid out as 2x2 tetramers, so this is the standard method.
#' 
#' @param input_file Path to the processed plate data CSV file
#' @param output_dir Directory to save the results (default: current directory)
#' @return Data frame with tetramer averages using 2x2 method
process_tetramer_averages <- function(input_file = "processed_plate_data.csv", output_dir = ".") {
  
  # Read the processed data
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }
  
  cat("Reading data from:", input_file, "\n")
  data <- read_csv(input_file, show_col_types = FALSE)
  
  # Calculate tetramer averages using 2x2 method (the standard layout)
  tetramer_data <- calculate_tetramer_averages_2x2(data)
  
  # Save results to specified directory
  output_file <- file.path(output_dir, "tetramer_averages_2x2.csv")
  write_csv(tetramer_data, output_file)
  
  cat("Tetramer averages saved to:", output_file, "\n")
  cat("Processed", nrow(tetramer_data), "tetramers\n")
  
  # Print summary statistics
  cat("\nSummary of tetramer averages:\n")
  print(summary(tetramer_data$readout_value))
  
  # Show distribution by gene
  gene_summary <- tetramer_data %>%
    group_by(gene) %>%
    summarise(
      count = n(),
      mean_readout = mean(readout_value, na.rm = TRUE),
      median_readout = median(readout_value, na.rm = TRUE),
      sd_readout = sd(readout_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_readout))
  
  cat("\nSummary by gene:\n")
  print(gene_summary)
  
  return(tetramer_data)
}

# Example usage:
# If running this script directly
if (!interactive()) {
  # Process using 2x2 method (standard layout)
  cat("Processing 2x2 tetramer averages...\n")
  tetramer_2x2 <- process_tetramer_averages()
  
  cat("\nTetramer averaging complete!\n")
}