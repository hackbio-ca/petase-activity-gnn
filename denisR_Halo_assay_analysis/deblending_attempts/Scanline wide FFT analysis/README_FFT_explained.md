# Enhanced FFT Explanation Script for HAIP Colony Analysis

## Overview

The enhanced `R/FFT_explained.R` script provides comprehensive visualization of how Fast Fourier Transform (FFT) improves colony analysis in the HAIP pipeline. It now features **dual-frequency detection** to identify both intra-tetramer and inter-tetramer spacings, **individual colony center mapping**, and **robust per-colony area calculation**.

## Key Enhancements

### 🔬 **Dual-Frequency Analysis**
- **Intra-tetramer frequency**: Detects spacing between colonies within a tetramer group
- **Inter-tetramer frequency**: Detects spacing between tetramer groups
- **Gitter-seeded detection**: Uses .dat file colony positions to seed frequency search bands
- **Harmonic redundancy filtering**: Prevents double-counting of harmonically related peaks

### 🎯 **Individual Colony Detection** 
- **Hilbert envelope analysis**: Uses envelope of narrowband-reconstructed signals to find colony centers
- **Local maxima detection**: Identifies individual colony positions rather than tetramer centers
- **Per-colony integration**: Calculates area for each colony using local minima bounds
- **Adaptive bounds**: Integration regions adapt to colony size and spacing

### 📈 **Robust Signal Processing**
- **Baseline detrending**: Running median or LOESS detrending removes low-frequency artifacts
- **Enhanced validation**: Automatic fallback to y-projection if gitter analysis fails
- **Power spectrum analysis**: Uses power spectrum instead of amplitude for better peak detection
- **Narrowband reconstruction**: Isolates specific frequency components for visualization

## Usage

```r
# Source the enhanced script
source("R/FFT_explained.R")

# Basic usage with enhanced dual-frequency detection
result <- fft_explained(
  washed_path = "path/to/washed_plate.jpg",
  background_path = "path/to/background_plate.jpg",
  gitter_dat_path = "path/to/plate.jpg.dat",
  out_png = "results/figs/fft_explained.png",
  out_csv = "results/fft_explained_metrics.csv"
)

# Results now include individual colony data
cat("Detected", length(result$colony_mapping$colony_centers), "individual colonies")
cat("Per-colony areas range:", range(result$per_colony_results$colony_areas$area))
```

## Enhanced Output Files

### 1. Visualization (PNG) - 3 Enhanced Panels

**Panel 1 (Top)**: Multi-signal comparison
- Raw washed signal (blue)
- Background signal (gray) 
- Background-adjusted signal (red)
- **NEW**: Detrended signal (dark red)
- **Individual colony centers** marked with green dashed lines (C1, C2, ...)
- Shaded integration area for selected colony with tight bounds

**Panel 2 (Middle)**: Frequency component analysis
- Windowed detrended signal (black)
- **NEW**: Intra-tetramer component overlay (red solid)
- **NEW**: Inter-tetramer component overlay (orange solid)
- **NEW**: Envelope traces (dotted lines) showing colony detection method
- Dynamic legend with actual frequency values

**Panel 3 (Bottom)**: Enhanced power spectrum
- **Power spectrum** (more robust than amplitude spectrum)
- **Dual frequency peaks** labeled:
  - `f_intra` (red) with period in pixels
  - `f_inter` (orange) with period in pixels  
- DC component (gray)
- Clear frequency-to-spacing mapping

### 2. Main Metrics CSV
Enhanced columns:
```
y_row, f_intra, period_intra, f_inter, period_inter, 
n_colonies_detected, colony_index_used, x_start, x_end, area, 
detrend_method, blur_sigma, notes
```

### 3. Per-Colony CSV (NEW)
Individual colony data:
```
colony_id, center_x, left_bound, right_bound, area
```

## Technical Implementation

### Enhanced Signal Processing Pipeline

1. **Load & Align**: Images loaded, gitter data parsed for row selection
2. **Scanline Extraction**: Horizontal line through first row center
3. **Background Adjustment**: `adjusted = washed - background + mean(background)`
4. **🆕 Robust Detrending**: Running median baseline removal
5. **FFT Preparation**: DC removal + Hann windowing
6. **🆕 Dual-Frequency Detection**: Gitter-seeded search in intra/inter bands
7. **🆕 Narrowband Reconstruction**: Component isolation + Hilbert envelopes
8. **🆕 Individual Colony Mapping**: Envelope maxima → colony centers
9. **🆕 Per-Colony Integration**: Local minima bounds for each colony

### Frequency Analysis Details

**Gitter-Seeded Detection**:
- Analyzes first row colony positions from .dat file
- Estimates intra-pair and inter-pair distances
- Converts to frequency seeds: `f = 1/distance`
- Searches ±25% around each seed frequency

**Peak Detection**:
- Uses power spectrum: `PSD = |FFT(signal)|²`
- Finds local maxima within search bands
- Filters harmonically redundant peaks
- Returns strongest non-harmonic peaks

**Component Reconstruction**:
- Keeps ±2 frequency bins around each peak
- Inverse FFT to time domain
- Hilbert transform for envelope calculation
- Local maxima in envelope = colony centers

### Per-Colony Area Calculation

**Adaptive Integration Bounds**:
- For each colony center, find nearest local minima left/right
- Integration bounds = [left_minimum, right_minimum]
- Trapezoidal integration over detrended signal
- Bounds adapt to colony size and spacing variations

**Validation & Fallback**:
- Validates intra frequency power vs background noise
- Auto-switches to y-projection if gitter analysis fails
- Logs all fallback decisions in CSV notes

## Dependencies

**Required R packages**:
- `imager`: Image loading and processing
- `ggplot2`: Enhanced visualization
- `gridExtra`: Multi-panel figure arrangement  
- `readr`: CSV file I/O
- `dplyr`: Data manipulation
- **🆕 `pracma`**: Peak finding and Hilbert transform

Install pracma: `install.packages("pracma")`

## Expected Results

### Console Output Example:
```
Enhanced FFT Colony Analysis Explanation
=======================================
Using gitter file to determine first row center...
Found 16 colonies in first row, mean y = 371.0
Estimated spacings - Intra: 45.2 px (f=0.0221), Inter: 98.7 px (f=0.0101)
Search bands - Intra: [0.0166, 0.0276], Inter: [0.0076, 0.0127]
Detected frequencies - Intra: 0.0218 (period=45.9), Inter: 0.0098 (period=102.1)
Found 14 individual colony centers using envelope detection
Calculated areas for 14 colonies
Area range: 1250.3 to 2847.9
```

### Validation Checks:
- **Two labeled peaks** in bottom spectrum (intra & inter)
- **Individual colony markers** (not tetramer markers) in top plot
- **Tight integration bounds** around selected colony
- **Per-colony CSV** with 14 rows (matching detected colonies)
- **Area sensitivity** to colony_index_for_area parameter

## Integration with HAIP Pipeline

**Backward Compatibility**: 
- Same function interface as original version
- Existing test scripts work with enhanced features
- Original output files still generated (with enhanced data)

**New Capabilities**:
- Colony-level intensity quantification for tetramer analysis
- Validation of FFT-based colony detection parameters
- Educational visualization of dual-frequency tetramer structure
- Per-colony area comparison across experimental conditions

The enhanced script provides the theoretical foundation and practical validation needed for robust FFT-based colony analysis in the HAIP pipeline.