# GMM Scanline Analysis Implementation

This implementation provides a complete GMM-based scanline colony quantification system as specified. The code consists of two main R scripts and supporting files.

## Files Created

### Core Implementation
- **`R/GMM_scanline.R`** - Main analysis functions
- **`R/GMM_explained.R`** - Visualization functions
- **`test_gmm_scanline.R`** - Test script for BHET25 dataset
- **`example_usage.R`** - Usage examples and documentation

### Output Directories
- **`results/scanline/`** - CSV output files
- **`results/figs/colony/`** - Plot output files

## Key Features

### GMM_scanline.R Functions

1. **`gmm_scanline()`** - Main analysis function
   - Reads washed, background, and gitter data
   - Crops images using gitter bounding boxes + margin
   - Extracts scanline and applies background correction
   - Detrends using running median
   - Defines colony windows from gitter centers
   - Fits 2-component GMM per window
   - Computes responsibility-weighted area
   - Applies presence/absence criteria
   - Saves CSV results

2. **Helper Functions:**
   - `read_gitter_dat()` - Parse gitter .dat files
   - `compute_global_crop()` - Calculate crop boundaries
   - `extract_scanline()` - Extract horizontal scanline
   - `detrend_runmed()` - Running median detrending
   - `define_windows_from_centers()` - Create analysis windows
   - `fit_gmm_window()` - Fit 2-component GMM
   - `compute_weighted_area()` - Calculate weighted area
   - `presence_check()` - Apply presence criteria

### GMM_explained.R Functions

1. **`gmm_explained()`** - Main visualization function
   - Creates overview plot of full scanline
   - Generates per-colony detailed plots (3 panels each)
   - Panel A: Full scanline context
   - Panel B: Window zoom with GMM fit
   - Panel C: Responsibility weights

2. **Helper Functions:**
   - `rebuild_scanline_for_plotting()` - Reconstruct analysis for plots
   - `plot_overview()` - Create scanline overview
   - `plot_colony()` - Create 3-panel colony plots
   - `reconstruct_gmm_for_plotting()` - Rebuild GMM for visualization

## Algorithm Implementation

### 1. Image Processing
- Loads JPG images using `imager::load.image()`
- Converts to grayscale using standard luminance formula
- Crops both images to union of gitter bounding boxes + margin

### 2. Scanline Extraction
- Uses mean y-coordinate of target gitter row
- Extracts horizontal line from cropped images
- Background correction: `adjusted = washed - background + mean(background)`

### 3. Detrending
- Running median with window size = 3.5 × median center spacing
- Clamped to [41, 201] pixels, forced to odd number
- Baseline removal: `detrended = adjusted - baseline`

### 4. Window Definition
- Sort gitter centers for target row
- Use midpoints between adjacent centers as initial bounds
- Snap bounds to local minima within ±snap_radius
- Apply distance constraints [r_intra_min, r_intra_max]
- Reject windows touching crop edges

### 5. GMM Fitting
- Uses `mixtools::normalmixEM()` with k=2 components
- Initializes means with 20th/80th percentiles
- Identifies colony component as higher mean
- Falls back to `mclust` if mixtools fails
- Computes posterior responsibilities γᵢ

### 6. Metrics Calculation
- **Area**: Σ γᵢ × (adjustedᵢ - baselineᵢ)
- **Lift**: Σ γᵢ × (adjustedᵢ - baselineᵢ) / Σ γᵢ
- **SNR**: μ_lift / MAD(adjusted_window - baseline_window)
- **Responsibility mean**: mean(γᵢ)

### 7. Presence Criteria
Colony present if ALL of:
- area ≥ area_min
- SNR ≥ snr_min  
- resp_mean ≥ resp_min
- Window doesn't touch crop edges

## Output Files

### gmm_scanline_per_colony.csv
Per-colony metrics:
- `colony_id`, `present`, `center_x`, `left`, `right`, `width_px`
- `area_weighted`, `mu_lift`, `snr`, `resp_mean`
- `gmm_mu_bg`, `gmm_sd_bg`, `gmm_pi_bg`
- `gmm_mu_col`, `gmm_sd_col`, `gmm_pi_col`
- `touches_edge`, `notes`, `y_row`, `row_index`, `source_centers`

### gmm_scanline_overview.csv
Run metadata:
- `y_row`, `row_index`, `n_colonies_total`, `n_colonies_present`
- Crop boundaries: `crop_left`, `crop_right`, `crop_top`, `crop_bottom`
- Parameters: `detrend_k`, `area_min`, `snr_min`, `resp_min`, etc.
- File paths

### Visualization Plots
- `scanline_overview.png` - Full scanline with present colonies
- `C1.png`, `C2.png`, ... - Individual colony plots (3 panels each)

## Default Parameters for BHET25

Optimized for row 1 (captures C9 & C10 tails):
```r
gmm_scanline(
  washed_path, background_path, gitter_dat_path,
  row_index = 1,
  margin = 25,
  snap_radius = 30,
  r_intra_min = 110,
  r_intra_max = 360,
  area_min = 2.0,
  snr_min = 0.6,
  resp_min = 0.25
)
```

## Required R Packages
- `imager` - Image loading and processing
- `mixtools` - Gaussian mixture model fitting
- `readr` - CSV file I/O
- `dplyr` - Data manipulation
- `ggplot2` - Plotting
- `gridExtra` - Multi-panel plots

## Usage

1. Load the functions:
```r
source("R/GMM_scanline.R")
source("R/GMM_explained.R")
```

2. Run analysis:
```r
result <- gmm_scanline(washed_path, background_path, gitter_dat_path)
```

3. Generate plots:
```r
viz_result <- gmm_explained(washed_path, background_path)
```

4. Check outputs in `results/scanline/` and `results/figs/colony/`

## Test Script

Run `test_gmm_scanline.R` to test with BHET25 dataset. The script will:
- Load test images and gitter data
- Run analysis with optimized parameters
- Generate all visualization plots
- Report number of present colonies found

This implementation follows the exact specifications provided and should work seamlessly with the existing HAIP codebase structure.