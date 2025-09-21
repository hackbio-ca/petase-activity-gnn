# Example usage of GMM scanline analysis
# 
# This script demonstrates how to use the GMM_scanline.R and GMM_explained.R
# functions for colony quantification.

# Load the functions
source("R/GMM_scanline.R")
source("R/GMM_explained.R")

# Example 1: Basic usage with BHET25 dataset
# Set up file paths
washed_path <- "Plate Images/washed/washed_BHET25_2d_6_1.JPG"
background_path <- "Plate Images/background/background_BHET25_6_1.JPG"
gitter_dat_path <- "Plate Images/unwashed/BHET25_2d_6_1.JPG.dat"

# Run scanline analysis for row 1 with default parameters optimized for BHET25
result <- gmm_scanline(
  washed_path = washed_path,
  background_path = background_path,
  gitter_dat_path = gitter_dat_path,
  out_dir = "results/scanline",
  row_index = 1,
  margin = 25,
  snap_radius = 30,
  r_intra_min = 110,
  r_intra_max = 360,
  area_min = 2.0,
  snr_min = 0.6,
  resp_min = 0.25
)

# Generate visualization plots
viz_result <- gmm_explained(
  washed_path = washed_path,
  background_path = background_path,
  scanline_csv = "results/scanline/gmm_scanline_per_colony.csv",
  overview_csv = "results/scanline/gmm_scanline_overview.csv",
  out_dir = "results/figs/colony"
)

# Example 2: Custom parameters for different scenarios
# If tails are clipped, try these adjusted parameters:
result_adjusted <- gmm_scanline(
  washed_path = washed_path,
  background_path = background_path,
  gitter_dat_path = gitter_dat_path,
  out_dir = "results/scanline",
  row_index = 1,
  margin = 25,
  snap_radius = 40,        # Increased for better tail capture
  r_intra_min = 110,
  r_intra_max = 400,       # Increased for wider search
  area_min = 2.0,
  snr_min = 0.4,          # Relaxed for borderline colonies
  resp_min = 0.25
)

# Example 3: Loop over multiple rows (future extension)
# for (row_i in 1:3) {
#   result_row <- gmm_scanline(
#     washed_path = washed_path,
#     background_path = background_path,
#     gitter_dat_path = gitter_dat_path,
#     out_dir = paste0("results/scanline/row_", row_i),
#     row_index = row_i
#   )
# }

cat("Analysis complete. Check results/scanline/ for CSV files and results/figs/colony/ for plots.\n")