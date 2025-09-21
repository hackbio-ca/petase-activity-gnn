# Main Plate Analysis Script for Halo Assay Analysis
# This script processes plate images and generates analysis data similar to plate_data_subset.csv

# Load required libraries
library(imager)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)

# Load HAIP package (assuming it's installed)
# If not installed, uncomment the following lines:
# devtools::install("./HAIP")
library(HAIP)

#' Process a single plate image and extract intensity data
#' 
#' @param image_path Path to the plate image file
#' @param dat_path Path to the corresponding .dat file with colony coordinates
#' @param background_path Path to the background image (optional)
#' @param plate_id Plate identifier (e.g., "HA00009")
#' @param plasmid_id Plasmid identifier (e.g., 152)
#' @return Data frame with colony analysis results
process_plate_image <- function(image_path, dat_path, background_path = NULL, plate_id = NULL, plasmid_id = NULL) {
  
  cat("Processing:", basename(image_path), "\n")
  
  # Load the main image
  img <- load.image(image_path)
  
  # Check if .dat file exists
  if (!file.exists(dat_path)) {
    cat("Warning: .dat file not found for", basename(image_path), "\n")  

    return(NULL)
  }
  
  # Load background image if provided
  bg_img <- NULL
  if (!is.null(background_path) && file.exists(background_path)) {
    bg_img <- load.image(background_path)
  }
  
  # Initialize colony data from .dat file
  colony_data <- initialize_colony_data(dat_path, img)
  
  # Create grid boxes for analysis
  colony_data <- create_grid_boxes(colony_data)
  
  # Calculate average intensity for each colony
  colony_data <- calculate_average_intensity(img, colony_data, prefix = "main")
  
  # If background image is available, calculate background intensity
  if (!is.null(bg_img)) {
    colony_data <- calculate_average_intensity(bg_img, colony_data, prefix = "bg")
    
    # Calculate edge intensity difference for normalization
    edge_result <- calculate_edge_intensity(img, bg_img)
    
    # Apply background correction
    colony_data$readout_value <- colony_data$main_mean_intensity - colony_data$bg_mean_intensity - edge_result$difference
  } else {
    # If no background, use raw intensity
    colony_data$readout_value <- colony_data$main_mean_intensity
  }
  
  # Handle NA values in readout_value with detailed logging
  na_indices <- which(is.na(colony_data$readout_value))
  na_count <- length(na_indices)
  
  if (na_count > 0) {
    cat("Warning: Found", na_count, "NA values in readout_value. Details:\n")
    
    # Create log entry with detailed information
    na_details <- data.frame(
      plate = plate_id,
      row = colony_data$row[na_indices],
      column = colony_data$col[na_indices],
      main_intensity = colony_data$main_mean_intensity[na_indices],
      bg_intensity = if("bg_mean_intensity" %in% names(colony_data)) colony_data$bg_mean_intensity[na_indices] else NA,
      timestamp = Sys.time()
    )
    
    # Print details to console
    for (i in 1:nrow(na_details)) {
      cat("  - Plate:", na_details$plate[i], "Row:", na_details$row[i], 
          "Column:", na_details$column[i], "\n")
    }
    
    # Write to log file
    log_file <- "na_values_log.csv"
    if (file.exists(log_file)) {
      existing_log <- read_csv(log_file, show_col_types = FALSE)
      combined_log <- bind_rows(existing_log, na_details)
    } else {
      combined_log <- na_details
    }
    write_csv(combined_log, log_file)
    
    # Replace NA values with 0
    colony_data$readout_value[na_indices] <- 0
    cat("  Replaced", na_count, "NA values with 0. Details saved to", log_file, "\n")
  }
  
  # Add metadata columns
  colony_data$plate <- plate_id
  colony_data$plasmid <- plasmid_id
  colony_data$normalization_method <- "raw"
  colony_data$colony_size <- colony_data$size
  colony_data$date_entered <- Sys.time()
  
  # Rename columns to match expected format
  colony_data$column <- colony_data$col
  
  # Select and arrange columns to match plate_data_subset.csv format
  result <- colony_data %>%
    select(plate, plasmid, column, row, normalization_method, readout_value, colony_size, date_entered) %>%
    arrange(row, column)
  
  return(result)
}

#' Get gene mapping for colonies based on position
#' 
#' This function maps colony positions to gene names using the reference plate data
#' 
#' @param data Data frame with plate, row and column positions
#' @return Data frame with gene column added
add_gene_mapping <- function(data) {
  # Load the reference plate data that contains the correct gene mappings
  plate_data_file <- "plate_data_subset.csv"
  
  if (!file.exists(plate_data_file)) {
    cat("Warning: Reference plate data file not found:", plate_data_file, "\n")
    cat("Using fallback mapping...\n")
    data$gene <- "Unknown"
    return(data)
  }
  
  # Read the reference data
  cat("Loading reference gene mapping from:", plate_data_file, "\n")
  reference_data <- read_csv(plate_data_file, show_col_types = FALSE)
  
  # Create a lookup table with unique plate/row/column to gene mappings
  gene_lookup <- reference_data %>%
    select(plate, row, column, gene) %>%
    distinct()
  
  # Join with the input data to get gene assignments
  data_with_genes <- data %>%
    left_join(gene_lookup, by = c("plate", "row", "column"))
  
  # Handle any remaining unmapped positions
  unmapped_count <- sum(is.na(data_with_genes$gene))
  if (unmapped_count > 0) {
    cat("Warning:", unmapped_count, "colonies could not be mapped to genes\n")
    data_with_genes$gene[is.na(data_with_genes$gene)] <- "Unknown"
  }
  
  return(data_with_genes)
}

#' Extract plate information from filename
#' 
#' @param filename Image filename
#' @return List with plate info
extract_plate_info <- function(filename) {
  # Extract information from filename like "BHET25_2d_6_1.JPG"
  # Map to existing plate IDs from the reference data
  
  base_name <- tools::file_path_sans_ext(basename(filename))
  parts <- str_split(base_name, "_")[[1]]
  
  # Create a mapping table for the available images to existing plate IDs
  # Based on available reference plates: HA00009, HA00021, HA00033, HA00045, HA00057, HA00069
  filename_to_plate <- list(
    "BHET25_2d_6_1" = "HA00009",
    "BHET25_2d_6_2" = "HA00021", 
    "BHET25_2d_6_3" = "HA00033",
    "BHET25_2d_7_4" = "HA00045",
    "BHET25_2d_7_5" = "HA00057",
    "BHET25_2d_7_6" = "HA00069"
  )
  
  plate_id <- filename_to_plate[[base_name]]
  
  if (is.null(plate_id)) {
    cat("Warning: No mapping found for", base_name, "using HA00009 as default\n")
    plate_id <- "HA00009"
  }
  
  return(list(
    plate_id = plate_id,
    plasmid_id = 152  # Default plasmid ID
  ))
}  

    "BHET25_2d_6_1" = "HA00009",  plate_images_dir <- file.path(base_dir, "Plate Images")

#' Main function to process all plate images
#' 
#' @param base_dir Base directory containing Plate Images folder
#' @return Combined data frame with all results
process_all_plates <- function(base_dir = ".") {
  
  plate_images_dir <- file.path(base_dir, "Plate Images")
  unwashed_dir <- file.path(plate_images_dir, "unwashed")
  washed_dir <- file.path(plate_images_dir, "washed")
  background_dir <- file.path(plate_images_dir, "background")
  
  # Get list of plate images (excluding gridded versions)
  unwashed_images <- list.files(unwashed_dir, pattern = "^BHET25.*\\.JPG$", full.names = TRUE)
  
  if (length(unwashed_images) == 0) {
    stop("No unwashed plate images found in: ", unwashed_dir)
  }
  
  cat("Found", length(unwashed_images), "plate images to process\n")
  
  all_results <- list()
  
  for (img_path in unwashed_images) {
    # Get corresponding .dat file
    dat_path <- paste0(img_path, ".dat")
    
    if (!file.exists(dat_path)) {
      cat("Warning: No .dat file found for", basename(img_path), "\n")
      next
    }
    
    # Get corresponding background image
    img_basename <- basename(img_path)
    # Remove "2d_" from the filename to match background naming pattern
    bg_basename <- paste0("background_", gsub("_2d_", "_", img_basename))
    bg_path <- file.path(background_dir, bg_basename)
    
    if (!file.exists(bg_path)) {
      cat("Warning: No background image found for", basename(img_path), "\n")
      bg_path <- NULL
    }
    
    # Extract plate information
    plate_info <- extract_plate_info(img_path)
    
    # Process the image
    result <- process_plate_image(
      image_path = img_path,
      dat_path = dat_path,
      background_path = bg_path,
      plate_id = plate_info$plate_id,
      plasmid_id = plate_info$plasmid_id
    )
    
    if (is.null(result) || nrow(result) == 0) {
      cat("Warning: No results for", basename(img_path), "\n")
      next
    }
    
    # Add gene mapping
    result <- add_gene_mapping(result)
    
    # Add unique ID (similar to the original data)
    result$id <- seq(from = nrow(result) * length(all_results) + 16129, 
                     length.out = nrow(result))
    
    all_results[[length(all_results) + 1]] <- result
  }
  
  if (length(all_results) == 0) {
    stop("No plates were successfully processed")
  }
  
  # Combine all results
  final_data <- bind_rows(all_results)
  
  # Reorder columns to match plate_data_subset.csv
  final_data <- final_data %>%
    select(id, gene, plate, plasmid, column, row, normalization_method, 
           readout_value, colony_size, date_entered)
  
  return(final_data)
}