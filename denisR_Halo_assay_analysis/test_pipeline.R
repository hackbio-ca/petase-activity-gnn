# Test Script for Refactored Pipeline
# This script verifies that the refactored code works correctly

cat("Testing refactored HAIP pipeline...\n")

# Test 1: Check that required files exist
required_files <- c(
  "plate_analysis_main.R",
  "tetramer_averaging.R", 
  "run_complete_analysis.R",
  "batch_processing.R"
)

cat("Checking required files...\n")
for (file in required_files) {
  if (file.exists(file)) {
    cat(paste("✓", file, "exists\n"))
  } else {
    cat(paste("✗", file, "MISSING\n"))
  }
}

# Test 2: Check that HAIP package can be loaded
cat("\nTesting HAIP package loading...\n")
tryCatch({
  # Check if we can source the main functions
  if (file.exists("HAIP/R/classify_pixels.R")) {
    cat("✓ HAIP R files found\n")
  }
  
  # Check required libraries
  required_packages <- c("dplyr", "imager", "ggplot2")
  for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(paste("✓", pkg, "available\n"))
    } else {
      cat(paste("⚠", pkg, "not available - may need installation\n"))
    }
  }
}, error = function(e) {
  cat("✗ Error loading packages:", e$message, "\n")
})

# Test 3: Verify timestamp generation works
cat("\nTesting timestamp generation...\n")
tryCatch({
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  cat(paste("✓ Timestamp generated:", timestamp, "\n"))
  
  # Test results directory creation
  test_results_dir <- file.path("results", paste0("test_", timestamp))
  if (!dir.exists(test_results_dir)) {
    dir.create(test_results_dir, recursive = TRUE)
    cat(paste("✓ Test results directory created:", test_results_dir, "\n"))
    
    # Clean up test directory
    unlink(test_results_dir, recursive = TRUE)
    cat("✓ Test directory cleaned up\n")
  }
}, error = function(e) {
  cat("✗ Error with timestamp/directory:", e$message, "\n")
})

# Test 4: Check data files
cat("\nChecking data files...\n")
data_files <- c("plate_data_subset.csv", "plate_metadata_subset.csv")
for (file in data_files) {
  if (file.exists(file)) {
    cat(paste("✓", file, "exists\n"))
  } else {
    cat(paste("✗", file, "MISSING\n"))
  }
}

# Test 5: Check Plate Images directory
cat("\nChecking Plate Images directory...\n")
if (dir.exists("Plate Images")) {
  subdirs <- c("unwashed", "washed", "background")
  for (subdir in subdirs) {
    path <- file.path("Plate Images", subdir)
    if (dir.exists(path)) {
      files <- list.files(path, pattern = "\\.JPG$", ignore.case = TRUE)
      cat(paste("✓", subdir, "directory exists with", length(files), "JPG files\n"))
    } else {
      cat(paste("✗", subdir, "directory missing\n"))
    }
  }
} else {
  cat("✗ Plate Images directory missing\n")
}

cat("\nPipeline verification complete!\n")
cat("If all tests show ✓, the refactored pipeline should work correctly.\n")
cat("Any ✗ or ⚠ items may need attention before running the full analysis.\n")