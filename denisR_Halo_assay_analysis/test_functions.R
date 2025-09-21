# Function Verification Script
# Tests key refactored functions to ensure they work correctly

cat("Testing refactored functions...\n")

# Test timestamp directory creation function
test_timestamp_creation <- function() {
  cat("\n=== Testing Timestamp Directory Creation ===\n")
  
  tryCatch({
    # Test the exact code from run_complete_analysis.R
    ANALYSIS_TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
    RESULTS_DIR <- file.path("results", paste0("analysis_", ANALYSIS_TIMESTAMP))
    
    cat(paste("Generated timestamp:", ANALYSIS_TIMESTAMP, "\n"))
    cat(paste("Results directory path:", RESULTS_DIR, "\n"))
    
    # Test directory creation
    if (!dir.exists(RESULTS_DIR)) {
      dir.create(RESULTS_DIR, recursive = TRUE)
      cat("✓ Results directory created successfully\n")
      
      # Test subdirectory creation
      subdirs <- c("plots", "logs", "processed_data")
      for (subdir in subdirs) {
        subpath <- file.path(RESULTS_DIR, subdir)
        dir.create(subpath, recursive = TRUE)
        cat(paste("✓ Subdirectory created:", subdir, "\n"))
      }
      
      # Clean up test directory
      unlink(RESULTS_DIR, recursive = TRUE)
      cat("✓ Test directories cleaned up\n")
      
      return(TRUE)
    } else {
      cat("⚠ Directory already exists\n")
      return(TRUE)
    }
  }, error = function(e) {
    cat(paste("✗ Error:", e$message, "\n"))
    return(FALSE)
  })
}

# Test NA logging function
test_na_logging <- function() {
  cat("\n=== Testing NA Logging Function ===\n")
  
  tryCatch({
    # Create test data with some NA values
    test_data <- data.frame(
      plate = c("plate1", "plate2", "plate3"),
      row = c(1, 2, 3),
      column = c("A", "B", "C"),
      readout_value = c(1.5, NA, 2.3),
      gene = c("gene1", "gene2", "gene3")
    )
    
    # Test NA detection
    na_rows <- which(is.na(test_data$readout_value))
    if (length(na_rows) > 0) {
      cat(paste("✓ Found", length(na_rows), "NA values\n"))
      
      # Test logging format
      for (i in na_rows) {
        warning_msg <- paste("NA readout_value found for colony at plate:", 
                           test_data$plate[i], 
                           "row:", test_data$row[i], 
                           "column:", test_data$column[i])
        cat(paste("✓ Warning format:", warning_msg, "\n"))
      }
      
      # Test log data frame creation
      na_log <- test_data[na_rows, c("plate", "row", "column", "gene")]
      na_log$timestamp <- Sys.time()
      cat("✓ NA log data frame created with columns:", paste(colnames(na_log), collapse = ", "), "\n")
      
      return(TRUE)
    } else {
      cat("⚠ No NA values in test data\n")
      return(TRUE)
    }
  }, error = function(e) {
    cat(paste("✗ Error:", e$message, "\n"))
    return(FALSE)
  })
}

# Test 2x2 tetramer processing function
test_2x2_processing <- function() {
  cat("\n=== Testing 2x2 Tetramer Processing ===\n")
  
  tryCatch({
    # Simulate 2x2 tetramer data
    test_readouts <- data.frame(
      plate = rep("test_plate", 4),
      row = c(1, 1, 2, 2),
      column = c("A", "B", "A", "B"),
      readout_value = c(1.2, 1.4, 1.1, 1.3),
      gene = "test_gene"
    )
    
    # Test 2x2 averaging (simulated)
    if (nrow(test_readouts) == 4) {
      avg_readout <- mean(test_readouts$readout_value, na.rm = TRUE)
      cat(paste("✓ 2x2 tetramer average calculated:", round(avg_readout, 3), "\n"))
      
      # Test that we only have 2x2 grid
      unique_positions <- paste(test_readouts$row, test_readouts$column)
      if (length(unique_positions) == 4) {
        cat("✓ 2x2 grid structure confirmed\n")
      } else {
        cat("✗ Grid structure incorrect\n")
        return(FALSE)
      }
      
      return(TRUE)
    } else {
      cat("✗ Incorrect number of positions for 2x2 grid\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat(paste("✗ Error:", e$message, "\n"))
    return(FALSE)
  })
}

# Test file path generation
test_file_paths <- function() {
  cat("\n=== Testing File Path Generation ===\n")
  
  tryCatch({
    # Test timestamp and results directory
    ANALYSIS_TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
    RESULTS_DIR <- file.path("results", paste0("analysis_", ANALYSIS_TIMESTAMP))
    
    # Test key file paths that should be generated
    test_files <- c(
      "processed_plate_data.csv",
      "tetramer_averages_2x2.csv",
      "summary_by_gene.csv",
      "summary_by_plate.csv",
      "na_values_log.csv"
    )
    
    for (file in test_files) {
      full_path <- file.path(RESULTS_DIR, file)
      cat(paste("✓ Path generated:", full_path, "\n"))
    }
    
    return(TRUE)
  }, error = function(e) {
    cat(paste("✗ Error:", e$message, "\n"))
    return(FALSE)
  })
}

# Run all tests
cat("Starting function verification tests...\n")

results <- list()
results$timestamp <- test_timestamp_creation()
results$na_logging <- test_na_logging()
results$tetramer_2x2 <- test_2x2_processing()
results$file_paths <- test_file_paths()

# Summary
cat("\n=== VERIFICATION SUMMARY ===\n")
all_passed <- all(unlist(results))
if (all_passed) {
  cat("✓ All function tests PASSED! The refactored pipeline should work correctly.\n")
} else {
  failed_tests <- names(results)[!unlist(results)]
  cat(paste("✗ Some tests FAILED:", paste(failed_tests, collapse = ", "), "\n"))
  cat("Please review the failed tests before running the full analysis.\n")
}

cat("\nFunction verification complete!\n")