# Halo Assay Analysis Pipeline

This R-based analysis pipeline processes plate images from halo assays and generates comprehensive analysis results including tetramer averages and comparisons with existing data.

## Overview

The pipeline consists of several interconnected R scripts that:
1. Process plate images using the HAIP package
2. Extract colony intensity and size data
3. Calculate tetramer averages (bottom two values)
4. Generate comprehensive statistics and visualizations
5. Compare results with existing data

## Files Description

### Main Scripts
- `run_complete_analysis.R` - Master script that runs the entire pipeline
- `plate_analysis_main.R` - Core image processing and analysis functions
- `tetramer_averaging.R` - Calculates tetramer averages from processed data
- `batch_processing.R` - Orchestrates batch processing and generates reports
- `comparison_analysis.R` - Compares processed results with existing data

### Data Files (Required)
- `plate_data_subset.csv` - Existing plate data for comparison
- `plate_metadata_subset.csv` - Metadata for the plates
- `Plate Images/` directory containing:
  - `unwashed/` - Raw plate images and .dat files
  - `washed/` - Washed plate images (optional)
  - `background/` - Background images for normalization

## Prerequisites

### R Packages Required
- imager (for image processing)
- dplyr (for data manipulation)
- readr (for file I/O)
- ggplot2 (for visualization)
- gridExtra (for plot arrangements)
- stringr (for string manipulation)
- lubridate (for date handling)

### HAIP Package
The custom HAIP package must be available in the `./HAIP/` directory. The master script will attempt to install it automatically.

## Usage

### Quick Start (Recommended)
```r
# Run the complete analysis pipeline
source("run_complete_analysis.R")
```

This will:
1. Check and install required packages
2. Process all plate images
3. Calculate tetramer averages
4. Generate comprehensive reports and visualizations
5. Compare results with existing data

### Individual Components

#### 1. Process Plate Images Only
```r
source("plate_analysis_main.R")
results <- process_all_plates()
```

#### 2. Calculate Tetramer Averages Only
```r
source("tetramer_averaging.R")
tetramer_data <- process_tetramer_averages("processed_plate_data.csv", method = "1x4")
```

#### 3. Run Comparison Analysis Only
```r
source("comparison_analysis.R")
comparison <- compare_analysis_results()
```

## Tetramer Averaging Methods

The pipeline supports two tetramer grouping methods:

1. **1x4 Method**: Groups 4 consecutive columns in each row
2. **2x2 Method**: Groups colonies in 2x2 grid patterns

For each tetramer group, the pipeline:
- Identifies all colonies in the group
- Selects the two colonies with the lowest readout values
- Calculates their average as the representative value

## Output Files

### Data Files
- `processed_plate_data.csv` - Raw colony analysis results
- `tetramer_averages_1x4.csv` - Tetramer averages using 1x4 method
- `tetramer_averages_2x2.csv` - Tetramer averages using 2x2 method
- `summary_by_gene.csv` - Gene-level summary statistics
- `summary_by_plate.csv` - Plate-level summary statistics

### Reports
- `analysis_report.md` - Comprehensive analysis summary
- `comparison_report.md` - Detailed comparison with existing data

### Visualizations
- `plot_readout_distribution.png` - Distribution of readout values
- `plot_readout_by_gene.png` - Readout values by gene
- `plot_size_vs_readout.png` - Colony size vs readout correlation
- `plot_readout_by_plate.png` - Plate-to-plate comparison
- `plot_tetramer_by_gene.png` - Tetramer results by gene
- `comparison_distributions.png` - Distribution comparison
- `comparison_genes.png` - Gene-level comparison
- `comparison_correlation.png` - Correlation analysis
- `comparison_complete.png` - Three-way comparison

## Data Format

The pipeline expects:
- Plate images in JPG format
- Corresponding .dat files with colony coordinates and metadata
- Background images for normalization (optional)

### .dat File Format
Tab-separated files with columns:
- row, col, size, circularity, flags, x, y, xl, xr, yt, yb

### Output Data Format
Similar to `plate_data_subset.csv` with columns:
- id, gene, plate, plasmid, column, row, normalization_method, readout_value, colony_size, date_entered

## Troubleshooting

### Common Issues

1. **HAIP Package Not Found**
   - Ensure the `HAIP/` directory exists in the working directory
   - Install devtools: `install.packages("devtools")`

2. **Missing Image Files**
   - Check that `Plate Images/unwashed/` contains JPG files and corresponding .dat files
   - Verify file naming convention matches the expected pattern

3. **Memory Issues with Large Images**
   - Process images in smaller batches
   - Reduce image resolution if necessary

4. **Gene Mapping Errors**
   - Update the `add_gene_mapping()` function in `plate_analysis_main.R` to match your plate layout

### Error Messages

- **"Plate Images directory not found"**: Ensure the directory structure is correct
- **"No colonies processed"**: Check that .dat files exist and are properly formatted
- **"HAIP package required"**: Install the HAIP package or place it in the correct directory

## Customization

### Modifying Gene Mapping
Edit the `add_gene_mapping()` function in `plate_analysis_main.R` to match your specific plate layout.

### Adjusting Tetramer Grouping
Modify the grouping logic in `tetramer_averaging.R` to change how tetramers are defined.

### Custom Analysis Parameters
Adjust parameters in the processing functions:
- Edge width for background correction
- Statistical measures calculated
- Visualization parameters

## Support

For issues or questions:
1. Check the generated log messages for error details
2. Verify all prerequisite files and packages are available
3. Review the troubleshooting section above
4. Check individual script documentation for specific functions

## Version Information

This pipeline is designed for:
- R version 4.0 or higher
- HAIP package (custom image analysis)
- Standard CRAN packages as listed in prerequisites