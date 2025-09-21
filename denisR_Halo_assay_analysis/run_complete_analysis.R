# Master Analysis Script for Halo Assay - Enhanced with Deblending
# This script orchestrates the complete analysis pipeline with optional halo deblending
#
# NEW FEATURES:
# - Watershed-based halo deblending for overlapping colonies (especially tetramers)
# - Quality control metrics and validation
# - Method comparison between traditional and deblending approaches
# - Configurable parameters for different halo characteristics
# - Enhanced reporting with deblending-specific metrics
#
# USAGE:
# 1. Adjust DEBLENDING_PARAMS above based on your halo appearance
# 2. Set ENABLE_DEBLENDING = TRUE to use watershed deblending
# 3. Set COMPARE_METHODS = TRUE to compare with traditional analysis
# 4. Run the script as usual
#
# OUTPUT FILES:
# - processed_plate_data.csv: Main results (compatible with existing workflow)
# - processed_plate_data_deblended.csv: Detailed deblending metrics
# - deblending_qc_summary.csv: Quality control summary
# - processed_plate_data_traditional.csv: Traditional method results (if comparing)
#
# For more information on deblending, see: HAIP/inst/HALO_DEBLENDING.md

# Set working directory and suppress warnings for cleaner output
options(warn = -1)

# Create timestamped results directory
ANALYSIS_TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
RESULTS_DIR <- file.path("results", paste0("analysis_", ANALYSIS_TIMESTAMP))

# Create the results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}
dir.create(RESULTS_DIR, recursive = TRUE)

cat("==============================================\n")
cat("    HALO ASSAY ANALYSIS PIPELINE\n")
cat("     Enhanced with Deblending Support\n")
cat("==============================================\n")
cat("Starting complete analysis at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results will be saved to:", RESULTS_DIR, "\n\n")

# ===============================================================================
# CONFIGURATION SECTION - Adjust these settings as needed
# ===============================================================================

# Deblending Configuration
ENABLE_DEBLENDING <- TRUE  # Set to FALSE to use traditional analysis only

# Deblending Parameters (adjust based on your halo characteristics)
DEBLENDING_PARAMS <- list(
  polarity = "dark_on_light",    # "dark_on_light" or "light_on_dark"
  blur_sigma = 2,                # Denoising strength (1-4, higher = more blur)
  mask_quantile = 0.90,          # Halo threshold (0.85-0.95 for bright halos)
  max_radius_px = 120,           # Maximum halo radius in pixels
  hmin = 2,                      # Watershed suppression (1-5, higher = less detail)
  use_dog = TRUE,                # Use Difference of Gaussians filter
  snap_seeds = TRUE              # Snap seeds to distance transform maxima
)

# For tetramer plates with overlapping halos, consider these alternative settings:
# DEBLENDING_PARAMS$max_radius_px <- 100  # Smaller radius for tight tetramers
# DEBLENDING_PARAMS$mask_quantile <- 0.85  # More sensitive detection
# DEBLENDING_PARAMS$hmin <- 3              # More suppression of noise

# Analysis Options
COMPARE_METHODS <- TRUE        # Compare deblending vs traditional methods
GENERATE_QC_PLOTS <- TRUE      # Generate quality control visualizations
SAVE_SEGMENTATIONS <- FALSE    # Save segmentation overlays (large files)

# ===============================================================================

# Check for required libraries and install if needed
required_packages <- c("imager", "dplyr", "readr", "ggplot2", "gridExtra", "corrplot", "stringr", "lubridate", "tibble")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cran.r-project.org/")
}

# Load all required libraries
suppressPackageStartupMessages({
  library(imager)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(gridExtra)
  library(stringr)
  library(lubridate)
  library(tibble)
})

# Check if HAIP package is available and verify deblending functionality
if (!requireNamespace("HAIP", quietly = TRUE)) {
  cat("Installing HAIP package from local directory...\n")
  if (dir.exists("./HAIP")) {
    tryCatch({
      devtools::install("./HAIP")
      library(HAIP)
      cat("✓ HAIP package installed successfully\n")
    }, error = function(e) {
      cat("✗ Error installing HAIP package:", e$message, "\n")
      cat("Please install the HAIP package manually before running this script\n")
      stop("HAIP package required")
    })
  } else {
    cat("✗ HAIP directory not found. Please ensure HAIP package is available.\n")
    stop("HAIP package directory not found")
  }
} else {
  library(HAIP)
  cat("✓ HAIP package loaded\n")

  # Check if deblending functions are available
  if (ENABLE_DEBLENDING) {
    if (exists("deblend_colonies", mode = "function")) {
      cat("✓ Halo deblending functionality detected\n")
    } else {
      cat("⚠ Deblending functions not found - using traditional analysis only\n")
      ENABLE_DEBLENDING <- FALSE
    }
  }
}

# Verify required data files exist
required_files <- c("plate_data_subset.csv", "plate_metadata_subset.csv")
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  cat("✗ Missing required files:", paste(missing_files, collapse = ", "), "\n")
  stop("Required data files not found")
}

# Verify plate images directory exists
if (!dir.exists("Plate Images")) {
  cat("✗ Plate Images directory not found\n")
  stop("Plate Images directory required")
}

cat("✓ All prerequisites verified\n")

if (ENABLE_DEBLENDING) {
  cat("✓ Deblending mode enabled with parameters:\n")
  cat("  - Polarity:", DEBLENDING_PARAMS$polarity, "\n")
  cat("  - Blur sigma:", DEBLENDING_PARAMS$blur_sigma, "\n")
  cat("  - Mask quantile:", DEBLENDING_PARAMS$mask_quantile, "\n")
  cat("  - Max radius:", DEBLENDING_PARAMS$max_radius_px, "px\n")
} else {
  cat("ℹ Using traditional bounding box analysis\n")
}
cat("\n")

#' Enhanced plate processing with optional deblending
#'
#' @param base_dir Base directory for analysis
#' @param use_deblending Whether to use watershed deblending
#' @param deblend_params List of deblending parameters
#' @return Enhanced results with deblending metrics if enabled
process_all_plates_enhanced <- function(base_dir = ".", use_deblending = FALSE, deblend_params = NULL) {

  cat("Enhanced plate processing with deblending support...\n")

  # Source the original plate analysis script
  source("plate_analysis_main.R")

  if (!use_deblending) {
    # Use original processing
    return(process_all_plates(base_dir))
  }

  # Enhanced processing with deblending
  plate_images_dir <- file.path(base_dir, "Plate Images")
  unwashed_dir <- file.path(plate_images_dir, "unwashed")
  washed_dir <- file.path(plate_images_dir, "washed")
  background_dir <- file.path(plate_images_dir, "background")

  # Get list of plate images
  washed_images <- list.files(washed_dir, pattern = "^BHET25.*\\.JPG$", full.names = TRUE)

  all_results <- list()
  all_deblend_results <- list()
  qc_summary <- list()

  for (img_path in washed_images) {
    cat("Processing with deblending:", basename(img_path), "\n")

    # Get corresponding .dat file from unwashed directory
    img_basename <- basename(img_path)
    unwashed_path <- file.path(unwashed_dir, img_basename)
    dat_path <- paste0(unwashed_path, ".dat")

    if (!file.exists(dat_path)) {
      cat("Warning: No .dat file found for", basename(img_path), "\n")
      next
    }

    tryCatch({
      # Load images
      washed_img <- load.image(img_path)

      # Initialize colony data
      colony_data <- initialize_colony_data(dat_path)

      # Create grid boxes
      colony_data <- create_grid_boxes(colony_data)

      # Add centroid coordinates for deblending (extract from bounding boxes)
      colony_data$x <- (colony_data$new_xl + colony_data$new_xr) / 2
      colony_data$y <- (colony_data$new_yt + colony_data$new_yb) / 2
      colony_data$id <- seq_len(nrow(colony_data))

      # Calculate background intensity for normalization
      colony_data_with_bg <- calculate_average_intensity(washed_img, colony_data, prefix = "bg")

      # Calculate edge intensity for background correction
      edge_result <- tryCatch({
        # Find corresponding background image
        bg_basename <- paste0("background_", gsub("_2d_", "_", img_basename))
        bg_path <- file.path(background_dir, bg_basename)

        if (file.exists(bg_path)) {
          bg_img <- load.image(bg_path)
          calculate_edge_intensity(washed_img, bg_img)
        } else {
          list(difference = 0)
        }
      }, error = function(e) {
        list(difference = 0)
      })

      # Run deblending analysis
      deblend_results <- classify_pixels(
        img = washed_img,
        colony_data_with_bg = colony_data_with_bg,
        background_normalization = edge_result$difference,
        deblend = TRUE,
        polarity = deblend_params$polarity,
        blur_sigma = deblend_params$blur_sigma,
        mask_quantile = deblend_params$mask_quantile,
        max_radius_px = deblend_params$max_radius_px
      )

      # Extract results
      if (!is.null(deblend_results$deblend_metrics)) {
        metrics <- deblend_results$deblend_metrics

        # Convert to traditional format
        plate_info <- extract_plate_info(img_path)

        result <- data.frame(
          id = metrics$id,
          gene = rep("Unknown", nrow(metrics)),  # Will be mapped later
          plate = plate_info$plate_id,
          plasmid = plate_info$plasmid_id,
          column = colony_data$col[metrics$id],
          row = colony_data$row[metrics$id],
          normalization_method = "deblended",
          readout_value = metrics$mean_int,
          colony_size = metrics$npx,
          date_entered = Sys.time(),
          stringsAsFactors = FALSE
        )

        # Add gene mapping
        result <- add_gene_mapping(result)

        # Store deblending-specific results
        deblend_summary <- data.frame(
          plate = plate_info$plate_id,
          image = basename(img_path),
          total_colonies = nrow(metrics),
          successful_segmentations = sum(!is.na(metrics$mean_int)),
          mean_intensity = mean(metrics$mean_int, na.rm = TRUE),
          mean_pixel_count = mean(metrics$npx, na.rm = TRUE),
          leaky_colonies = sum(metrics$flag_leak > 0.1, na.rm = TRUE),
          stringsAsFactors = FALSE
        )

        all_results[[length(all_results) + 1]] <- result
        all_deblend_results[[length(all_deblend_results) + 1]] <- metrics
        qc_summary[[length(qc_summary) + 1]] <- deblend_summary

        cat("  ✓ Deblended", nrow(metrics), "colonies,",
            sum(!is.na(metrics$mean_int)), "successful\n")
      }

    }, error = function(e) {
      cat("  ✗ Error processing", basename(img_path), ":", e$message, "\n")
    })
  }

  if (length(all_results) == 0) {
    cat("No plates were successfully processed\n")
    return(data.frame())
  }

  # Combine results
  final_data <- bind_rows(all_results)
  final_deblend_data <- bind_rows(all_deblend_results)
  final_qc <- bind_rows(qc_summary)

  # Save deblending-specific outputs to timestamped results directory
  write_csv(final_deblend_data, file.path(RESULTS_DIR, "processed_plate_data_deblended.csv"))
  write_csv(final_qc, file.path(RESULTS_DIR, "deblending_qc_summary.csv"))

  cat("✓ Deblending results saved to", file.path(RESULTS_DIR, "processed_plate_data_deblended.csv"), "\n")
  cat("✓ QC summary saved to", file.path(RESULTS_DIR, "deblending_qc_summary.csv"), "\n")

  # Return data in traditional format for compatibility
  return(final_data)
}

# Run the complete analysis pipeline
tryCatch({

  # Step 1: Source all analysis scripts
  cat("STEP 1: Loading analysis scripts...\n")
  source("plate_analysis_main.R")
  source("tetramer_averaging.R")
  source("batch_processing.R")
  source("comparison_analysis.R")
  cat("✓ All scripts loaded successfully\n\n")

  # Step 2: Run plate image analysis
  cat("STEP 2: Processing plate images...\n")
  cat("This may take several minutes depending on the number of images...\n")

  if (ENABLE_DEBLENDING) {
    cat("Using enhanced deblending analysis...\n")
    plate_results <- process_all_plates_enhanced(".",
                                                 use_deblending = TRUE,
                                                 deblend_params = DEBLENDING_PARAMS)

    # Optional: Compare with traditional method
    if (COMPARE_METHODS) {
      cat("Running traditional analysis for comparison...\n")
      traditional_results <- process_all_plates_enhanced(".",
                                                        use_deblending = FALSE,
                                                        deblend_params = NULL)

      if (nrow(traditional_results) > 0) {
        write_csv(traditional_results, "processed_plate_data_traditional.csv")

        # Quick comparison metrics
        if (nrow(plate_results) > 0 && nrow(traditional_results) > 0) {
          deblend_mean <- mean(plate_results$readout_value, na.rm = TRUE)
          trad_mean <- mean(traditional_results$readout_value, na.rm = TRUE)

          cat("✓ Method comparison:\n")
          cat("  - Deblending mean intensity:", round(deblend_mean, 3), "\n")
          cat("  - Traditional mean intensity:", round(trad_mean, 3), "\n")
          cat("  - Difference:", round(abs(deblend_mean - trad_mean), 3), "\n")
        }
      }
    }
  } else {
    cat("Using traditional bounding box analysis...\n")
    plate_results <- process_all_plates(".")
  }

  if (nrow(plate_results) > 0) {
    cat("✓ Successfully processed", nrow(plate_results), "colonies from",
        length(unique(plate_results$plate)), "plates\n")

    # Save traditional format results to timestamped results directory
    write_csv(plate_results, file.path(RESULTS_DIR, "processed_plate_data.csv"))
    cat("✓ Results saved to", file.path(RESULTS_DIR, "processed_plate_data.csv"), "\n")
  } else {
    cat("✗ No colonies processed - check image files and .dat files\n")
    stop("Plate processing failed")
  }

  # Step 3: Calculate tetramer averages
  cat("\nSTEP 3: Calculating tetramer averages using 2x2 method...\n")

  # Process using 2x2 method (standard plate layout)
  tetramer_2x2 <- tryCatch({
    process_tetramer_averages(
      input_file = file.path(RESULTS_DIR, "processed_plate_data.csv"),
      output_dir = RESULTS_DIR
    )
  }, error = function(e) {
    cat("Warning: 2x2 tetramer method failed:", e$message, "\n")
    NULL
  })

  if (!is.null(tetramer_2x2)) {
    cat("✓ 2x2 tetramer averages calculated:", nrow(tetramer_2x2), "tetramers\n")
  }

  # Step 4: Run complete batch analysis
  cat("\nSTEP 4: Running comprehensive analysis...\n")
  analysis_results <- run_complete_analysis(".", tetramer_method = "2x2", results_dir = RESULTS_DIR)
  cat("✓ Comprehensive analysis completed\n")

  # Step 5: Run comparison analysis
  cat("\nSTEP 5: Comparing with existing data...\n")
  comparison_results <- compare_analysis_results(
    existing_data_file = "plate_data_subset.csv",
    processed_data_file = file.path(RESULTS_DIR, "processed_plate_data.csv"),
    tetramer_data_file = file.path(RESULTS_DIR, "tetramer_averages_2x2.csv"),
    results_dir = RESULTS_DIR
  )
  cat("✓ Comparison analysis completed\n")

  # Final summary
  cat("\n==============================================\n")
  cat("           ANALYSIS COMPLETE\n")
  cat("==============================================\n")
  cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  cat("FILES GENERATED IN", RESULTS_DIR, ":\n")
  cat("- processed_plate_data.csv: Raw colony analysis results\n")
  if (ENABLE_DEBLENDING && file.exists(file.path(RESULTS_DIR, "processed_plate_data_deblended.csv"))) {
    cat("- processed_plate_data_deblended.csv: Detailed deblending metrics\n")
    cat("- deblending_qc_summary.csv: Quality control summary for deblending\n")
  }
  if (file.exists(file.path(RESULTS_DIR, "tetramer_averages_2x2.csv"))) {
    cat("- tetramer_averages_2x2.csv: Tetramer averaged results (2x2 method)\n")
  }
  cat("- summary_by_gene.csv: Gene-level summary statistics\n")
  cat("- summary_by_plate.csv: Plate-level summary statistics\n")
  cat("- analysis_report.md: Comprehensive analysis report\n")
  cat("- comparison_report.md: Comparison analysis report\n")
  cat("- plot_*.png: Various visualization plots\n")
  cat("- comparison_*.png: Comparison visualization plots\n")
  if (file.exists(file.path(RESULTS_DIR, "na_values_log.csv"))) {
    cat("- na_values_log.csv: Log of all NA values found during processing\n")
  }
  cat("\n")

  cat("SUMMARY STATISTICS:\n")
  cat("- Total colonies processed:", nrow(plate_results), "\n")
  cat("- Unique plates:", length(unique(plate_results$plate)), "\n")
  cat("- Unique genes:", length(unique(plate_results$gene)), "\n")
  cat("- Mean readout value:", round(mean(plate_results$readout_value, na.rm = TRUE), 4), "\n")
  if (ENABLE_DEBLENDING) {
    cat("- Analysis method: Watershed deblending\n")
    if (file.exists("deblending_qc_summary.csv")) {
      qc_data <- read_csv("deblending_qc_summary.csv", show_col_types = FALSE)
      cat("- Deblending success rate:",
          round(mean(qc_data$successful_segmentations / qc_data$total_colonies) * 100, 1), "%\n")
      cat("- Average colonies per plate:", round(mean(qc_data$total_colonies), 1), "\n")
      cat("- Plates with leakage issues:", sum(qc_data$leaky_colonies > 0), "\n")
    }
  } else {
    cat("- Analysis method: Traditional bounding box\n")
  }
  if (!is.null(tetramer_2x2)) {
    cat("- Tetramers calculated (2x2):", nrow(tetramer_2x2), "\n")
  }

  cat("\nThe analysis has been completed successfully!\n")
  if (ENABLE_DEBLENDING) {
    cat("Deblending analysis provides improved accuracy for overlapping halos.\n")
    cat("Review deblending_qc_summary.csv for quality control metrics.\n")
  }
  cat("Check the generated files for detailed results.\n")

}, error = function(e) {
  cat("\n✗ ERROR in analysis pipeline:", e$message, "\n")
  cat("Check the error details above and ensure all prerequisites are met.\n")
  if (ENABLE_DEBLENDING) {
    cat("If deblending failed, try setting ENABLE_DEBLENDING <- FALSE and rerun.\n")
  }
  stop("Analysis pipeline failed")
})

#' Helper function to validate deblending results
#'
#' @param deblend_file Path to deblending results file
validate_deblending_results <- function(deblend_file = "processed_plate_data_deblended.csv") {
  if (!file.exists(deblend_file)) {
    cat("Deblending results file not found\n")
    return(FALSE)
  }

  tryCatch({
    data <- read_csv(deblend_file, show_col_types = FALSE)

    cat("\nDEBLENDING VALIDATION:\n")
    cat("- Colonies segmented:", nrow(data), "\n")
    cat("- Successful segmentations:", sum(!is.na(data$mean_int)), "\n")
    cat("- Success rate:", round(sum(!is.na(data$mean_int)) / nrow(data) * 100, 1), "%\n")
    cat("- Mean intensity range:", round(min(data$mean_int, na.rm = TRUE), 2),
        "to", round(max(data$mean_int, na.rm = TRUE), 2), "\n")
    cat("- Mean pixel count:", round(mean(data$npx, na.rm = TRUE), 0), "\n")

    # Check for potential issues
    high_leak <- sum(data$flag_leak > 0.1, na.rm = TRUE)
    if (high_leak > 0) {
      cat("⚠ Warning:", high_leak, "colonies show leakage flags\n")
    }

    return(TRUE)
  }, error = function(e) {
    cat("Error validating deblending results:", e$message, "\n")
    return(FALSE)
  })
}

# Run validation if deblending was used
if (ENABLE_DEBLENDING && exists("plate_results") && nrow(plate_results) > 0) {
  validate_deblending_results()
}

# Restore warnings
options(warn = 0)

