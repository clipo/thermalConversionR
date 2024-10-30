# Thermal Image Calibration Tools

A comprehensive R toolkit for calibrating and processing thermal camera images using various techniques.

## Features

* Multiple calibration methods:
  * Sequential Histogram Matching
  * CLAHE (Contrast Limited Adaptive Histogram Equalization)
  * MSR (Multi-Scale Retinex)
  * Bilateral Filtering
  * LCEP (Local Contrast Enhancement with Edge Preservation)

## Prerequisites

```R
install.packages(c("tiff", "EBImage", "grid"))
```

## Usage

### Basic Usage

```R
# Process all images using all methods
process_thermal_images("./images", "./processed_images")

# Use specific methods
process_thermal_images("./images", "./processed_images", 
                      methods = c("hist_matched", "clahe"))

# Use single method
process_thermal_images("./images", "./processed_images", 
                      methods = "bilateral")
```

### Available Methods

* `original`: Original thermal image (always included)
* `hist_matched`: Sequential histogram matching
* `clahe`: Contrast Limited Adaptive Histogram Equalization
* `msr`: Multi-Scale Retinex
* `bilateral`: Bilateral Filter
* `lcep`: Local Contrast Enhancement with Edge Preservation

## Output Structure

```
processed_images/
├── comparisons/
│   └── methods_comparison.png
└── processed/
    ├── hist_matched/
    ├── clahe/
    ├── msr/
    ├── bilateral/
    └── lcep/
```

* `comparisons/`: Contains side-by-side visual comparisons of all methods
* `processed/`: Contains processed images in method-specific subdirectories

## Method Details

### Sequential Histogram Matching
Matches histograms sequentially through the image sequence, using each processed frame as reference for the next. This helps maintain consistency across the thermal sequence.

### CLAHE
Enhances local contrast while limiting amplification to reduce noise. Particularly effective for thermal images with varying temperature distributions.

### Multi-Scale Retinex (MSR)
Enhances image details across multiple scales while maintaining temperature relationships. Uses scales of 2, 4, and 8 pixels.

### Bilateral Filter
Edge-preserving smoothing that reduces noise while maintaining sharp temperature boundaries. Uses spatial and range parameters for filtering.

### LCEP
Enhances local contrast while preserving edge information, useful for highlighting thermal features without introducing artifacts.

## Parameters

* `input_dir`: Directory containing TIFF format thermal images
* `output_dir`: Directory for processed outputs
* `methods`: Vector of method names or "all" for all methods

## Notes

* All methods preserve original temperature relationships
* Images are processed in alphanumeric order
* Original temperature ranges are maintained in saved TIFF files
* Visualization includes thermal colormap (black-purple-red-yellow-white)

## Error Handling

The script includes comprehensive error handling for:
* Missing input directory
* No TIFF files found
* Invalid method selection
* Processing errors for individual files

## Example

```R
# Load required libraries
library(tiff)
library(EBImage)
library(grid)

# Process thermal images
process_thermal_images(
  input_dir = "./thermal_sequence",
  output_dir = "./calibrated_output",
  methods = c("hist_matched", "clahe")
)