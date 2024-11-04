# Thermal Image Calibration Tools (v3.0)

A comprehensive R toolkit for calibrating and processing thermal camera images using various enhancement and calibration techniques.
![Example of image processing comparisons](methods_comparison.png)

## Version History
* v3.0 - Memory optimized CLAHE, enhanced smoothing, parallel processing support
* Previous versions focused on basic image processing and calibration methods

## Features

* Multiple calibration and enhancement methods:
  * Sequential Histogram Matching
  * CLAHE (Contrast Limited Adaptive Histogram Equalization) with memory optimization
  * MSR (Multi-Scale Retinex)
  * Bilateral Filtering
  * LCEP (Local Contrast Enhancement with Edge Preservation)
  * Adaptive Gamma Correction
  * Guided Filtering
  * Wavelet-based Fusion
  * Local Histogram Equalization
 
## Metadata Preservation with ExifTool
To ensure that the original metadata and EXIF headers are retained in the processed output files, this script integrates with `exiftool`.

### Installing ExifTool
- **Linux/Unix**: Install via your package manager (e.g., `sudo apt-get install exiftool`).
- **Windows**: Download and install ExifTool from the [official website](https://exiftool.org/).
- **macOS**: Use Homebrew to install: `brew install exiftool`.

### Ensuring ExifTool is Accessible
- Verify the installation by running:
  ```bash
  exiftool -ver

 If exiftool is not recognized, add it to your system’s PATH:
	- **Linux/Unix:** Add export PATH=$PATH:/path/to/exiftool to ~/.bashrc or ~/.zshrc.
	- **Windows:** Add the path to exiftool in Environment Variables under System Properties.

## Prerequisites

```R
# Install required packages
install.packages(c("tiff", "EBImage", "grid", "gridExtra", "parallel", "doParallel"))

# For EBImage (Bioconductor package)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage")
```

## Usage

### Basic Usage

```R
# Process all images using all methods
process_thermal_images("./images", "./processed_images")

# Use specific methods
process_thermal_images("./images", "./processed_images", 
                      methods = c("hist_matched", "clahe", "adaptive_gamma"))

# Use single method with parallel processing disabled
process_thermal_images("./images", "./processed_images", 
                      methods = "bilateral",
                      parallel = FALSE)
```

### Available Methods

* `original`: Original thermal image (always included)
* `hist_matched`: Sequential histogram matching
* `clahe`: Memory-optimized Contrast Limited Adaptive Histogram Equalization
* `msr`: Multi-Scale Retinex
* `bilateral`: Bilateral Filter
* `lcep`: Local Contrast Enhancement with Edge Preservation
* `adaptive_gamma`: Adaptive Gamma Correction
* `guided_filter`: Edge-Preserving Guided Filter
* `wavelet_fusion`: Multi-scale Wavelet Fusion
* `local_histogram`: Local Histogram Equalization

## Output Structure

```
processed_images/
├── comparisons/
│   ├── methods_comparison.png
│   └── processing_summary.txt
└── processed/
    ├── hist_matched/
    ├── clahe/
    ├── msr/
    ├── bilateral/
    ├── lcep/
    ├── adaptive_gamma/
    ├── guided_filter/
    ├── wavelet_fusion/
    └── local_histogram/
```

* `comparisons/`: Contains side-by-side visual comparisons and processing summary
* `processed/`: Contains processed images in method-specific subdirectories

## Performance Optimizations (v3.0)

* Memory-efficient CLAHE implementation
* Parallel processing support
* Optimized Gaussian filtering operations
* Improved tile processing with reduced memory footprint
* Smart garbage collection
* Limited visualization output (max 100 images) to prevent memory issues
* Progress tracking for long operations

## Method Details

### Memory-Optimized CLAHE (v3.0)
* Efficient tile-based processing with smooth blending
* Separated Gaussian filtering into 1D operations
* Adaptive tile size selection
* Multiple smoothing passes with memory management
* Contrast limiting with efficient redistribution

[Previous method descriptions remain the same...]

## Parameters

* `input_dir`: Directory containing TIFF format thermal images
* `output_dir`: Directory for processed outputs
* `methods`: Vector of method names or "all" for all methods
* `parallel`: Boolean to enable/disable parallel processing (default: TRUE)

## Memory Management (v3.0)

* Smart memory usage during large image processing
* Efficient matrix operations
* Garbage collection optimization
* Progress tracking with memory status
* Visualization limits to prevent memory overflow

## Notes

* All methods preserve original temperature relationships
* Images are processed in alphanumeric order
* Original temperature ranges are maintained in saved TIFF files
* Visualization includes thermal colormap (black-purple-red-yellow-white)
* Processing summary includes all images even if not all are visualized
* Memory-efficient processing suitable for large image sets

## Error Handling

The script includes comprehensive error handling for:
* Memory limitations
* Missing input directory
* No TIFF files found
* Invalid method selection
* Processing errors for individual files
* Matrix dimension mismatches

## Example

```R
# Load required libraries
library(tiff)
library(EBImage)
library(grid)
library(gridExtra)
library(parallel)
library(doParallel)

# Process thermal images with memory-optimized CLAHE
process_thermal_images(
  input_dir = "./thermal_sequence",
  output_dir = "./calibrated_output",
  methods = c("clahe", "adaptive_gamma"),
  parallel = TRUE
)
```
