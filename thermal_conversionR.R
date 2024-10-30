library(tiff)
library(EBImage)
library(OpenImageR)
library(grid)

process_thermal_images <- function(input_dir, output_dir, methods = "all") {
  # Available methods
  available_methods <- c("original", "hist_matched", "clahe", "msr", "bilateral", "lcep")
  
  # Process methods argument
  if (length(methods) == 1 && methods == "all") {
    selected_methods <- available_methods
  } else {
    # Validate methods
    invalid_methods <- setdiff(methods, available_methods)
    if (length(invalid_methods) > 0) {
      stop("Invalid methods: ", paste(invalid_methods, collapse = ", "), 
           "\nAvailable methods: ", paste(available_methods, collapse = ", "))
    }
    selected_methods <- methods
  }
  
  # Always include original if not selected
  if (!"original" %in% selected_methods) {
    selected_methods <- c("original", selected_methods)
  }
  
  # Check input directory
  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }
  
  # Get list of TIFF files
  tiff_files <- list.files(input_dir, 
                           pattern = "\\.tiff?$", 
                           full.names = TRUE, 
                           ignore.case = TRUE)
  
  # Check if any TIFF files were found
  if (length(tiff_files) == 0) {
    stop("No TIFF files found in directory: ", input_dir)
  }
  
  # Sort files to ensure consistent processing order
  tiff_files <- sort(tiff_files)
  cat("Found", length(tiff_files), "TIFF files\n")
  
  # Create output directories
  subdirs <- c("comparisons", "processed")
  sapply(file.path(output_dir, subdirs), 
         dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # Create method-specific subdirectories
  method_dirs <- setdiff(selected_methods, "original")
  sapply(file.path(output_dir, "processed", method_dirs), 
         dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate global min/max with proper handling of signed integers
  all_mins <- all_maxs <- numeric(length(tiff_files))
  for (i in seq_along(tiff_files)) {
    img <- readTIFF(tiff_files[i], all = TRUE, native = TRUE)[[1]]
    all_mins[i] <- min(img, na.rm = TRUE)
    all_maxs[i] <- max(img, na.rm = TRUE)
  }
  global_min <- min(all_mins)
  global_max <- max(all_maxs)
  
  # Helper function for thermal colorization
  colorize_thermal <- function(img) {
    normalized <- (img - global_min) / (global_max - global_min)
    normalized[normalized < 0] <- 0
    normalized[normalized > 1] <- 1
    colors <- colorRampPalette(c("black", "purple", "red", "yellow", "white"))(256)
    matrix(colors[floor(normalized * 255) + 1], nrow = nrow(img))
  }
  
  match_percentile_based <- function(source_img, reference_img) {
    # Normalize both images to [0,1] range
    src_norm <- (source_img - min(source_img)) / (diff(range(source_img)))
    ref_norm <- (reference_img - min(reference_img)) / (diff(range(reference_img)))
    
    # Calculate percentiles on normalized images
    ref_percentiles <- quantile(ref_norm, probs = seq(0, 1, 0.01))
    src_percentiles <- quantile(src_norm, probs = seq(0, 1, 0.01))
    
    # Map normalized source to normalized reference
    result_norm <- src_norm
    for(i in 1:length(src_norm)) {
      bin <- findInterval(src_norm[i], src_percentiles)
      bin <- min(max(bin, 1), length(ref_percentiles))
      result_norm[i] <- ref_percentiles[bin]
    }
    
    # Scale back to source range
    result <- result_norm * diff(range(source_img)) + min(source_img)
    
    return(result)
  }
  
  # Helper functions for other methods
  apply_msr <- function(img, scales = c(2, 4, 8)) {
    result <- matrix(0, nrow = nrow(img), ncol = ncol(img))
    img_norm <- (img - min(img)) / (max(img) - min(img))
    img_eb <- Image(img_norm)
    
    for (scale in scales) {
      gaussian <- gblur(img_eb, sigma = scale)
      result <- result + log(img_norm + 1) - log(as.matrix(gaussian) + 1)
    }
    result <- (result / length(scales))
    result <- (result - min(result)) / (max(result) - min(result)) * 
      (global_max - global_min) + global_min
    return(result)
  }
  
  apply_bilateral <- function(img, sigma_s = 2, sigma_r = 0.1) {
    img_norm <- (img - min(img)) / (max(img) - min(img))
    result <- img_norm
    half_size <- ceiling(3 * sigma_s)
    window_size <- 2 * half_size + 1
    
    spatial_kernel <- matrix(0, window_size, window_size)
    for(i in 1:window_size) {
      for(j in 1:window_size) {
        x <- i - half_size - 1
        y <- j - half_size - 1
        spatial_kernel[i,j] <- exp(-(x^2 + y^2)/(2 * sigma_s^2))
      }
    }
    spatial_kernel <- spatial_kernel / sum(spatial_kernel)
    
    padded <- matrix(0, nrow = nrow(img_norm) + 2*half_size, 
                     ncol = ncol(img_norm) + 2*half_size)
    padded[(half_size+1):(half_size+nrow(img_norm)), 
           (half_size+1):(half_size+ncol(img_norm))] <- img_norm
    
    for(i in 1:nrow(img_norm)) {
      if(i %% 10 == 0) cat("Processing row", i, "of", nrow(img_norm), "\n")
      for(j in 1:ncol(img_norm)) {
        window <- padded[i:(i+2*half_size), j:(j+2*half_size)]
        intensity_diff <- (window - img_norm[i,j])^2
        range_kernel <- exp(-intensity_diff/(2 * sigma_r^2))
        weights <- spatial_kernel * range_kernel
        result[i,j] <- sum(window * weights) / sum(weights)
      }
    }
    
    result <- result * (max(img) - min(img)) + min(img)
    return(result)
  }
  
  apply_lcep <- function(img, window_size = 7, alpha = 0.5) {
    img_norm <- (img - min(img)) / (max(img) - min(img))
    result <- matrix(0, nrow = nrow(img), ncol = ncol(img))
    
    for (i in 1:nrow(img)) {
      for (j in 1:ncol(img)) {
        r1 <- max(1, i - window_size%/%2)
        r2 <- min(nrow(img), i + window_size%/%2)
        c1 <- max(1, j - window_size%/%2)
        c2 <- min(ncol(img), j + window_size%/%2)
        
        window <- img_norm[r1:r2, c1:c2]
        local_mean <- mean(window)
        local_std <- sd(window)
        
        if (local_std > 0) {
          result[i,j] <- local_mean + alpha * (img_norm[i,j] - local_mean)
        } else {
          result[i,j] <- img_norm[i,j]
        }
      }
    }
    result <- result * (global_max - global_min) + global_min
    return(result)
  }
  
  # Process each image with all methods
  all_results <- list()
  reference_img <- NULL
  ref_mean <- NULL
  ref_std <- NULL
  
  # Track median values
  median_values <- data.frame(
    frame = integer(),
    original = numeric(),
    matched = numeric()
  )
  
  for (i in seq_along(tiff_files)) {
    img_path <- tiff_files[i]
    base_name <- tools::file_path_sans_ext(basename(img_path))
    cat("Processing:", base_name, "\n")
    
    tryCatch({
      curr_img <- readTIFF(img_path, all = TRUE, native = TRUE)[[1]]
      if (any(is.na(curr_img))) curr_img[is.na(curr_img)] <- global_min
      
      results <- list()
      results$original <- curr_img
      
      # Histogram Matching using percentile-based approach
      if ("hist_matched" %in% selected_methods) {
        if (i == 1) {
          # Store first image as reference
          reference_img <- curr_img
          results$hist_matched <- curr_img
        } else {
          # Match to first image using percentile matching
          results$hist_matched <- match_percentile_based(curr_img, reference_img)
        }
        
        # Track median values for this frame
        median_values <- rbind(median_values, 
                               data.frame(
                                 frame = i,
                                 original = median(curr_img),
                                 matched = median(results$hist_matched)
                               ))
        
        # Save histogram matched image
        writeTIFF(matrix(results$hist_matched, nrow=nrow(curr_img)),
                  file.path(output_dir, "processed", "hist_matched", 
                            paste0(base_name, "_hist_matched.tiff")),
                  bits.per.sample = 16,
                  compression = "none",
                  native = TRUE)
      }
      
      # CLAHE
      if ("clahe" %in% selected_methods) {
        scaled_img <- (curr_img - global_min) / (global_max - global_min)
        img_eb <- Image(scaled_img)
        clahe_img <- clahe(img_eb)
        results$clahe <- matrix(as.array(clahe_img) * (global_max - global_min) + global_min,
                                nrow = nrow(curr_img), 
                                ncol = ncol(curr_img))
        writeTIFF(results$clahe, 
                  file.path(output_dir, "processed", "clahe", 
                            paste0(base_name, "_clahe.tiff")),
                  bits.per.sample = 16,
                  compression = "none",
                  native = TRUE)
      }
      
      # MSR
      if ("msr" %in% selected_methods) {
        results$msr <- matrix(apply_msr(curr_img),
                              nrow = nrow(curr_img), 
                              ncol = ncol(curr_img))
        writeTIFF(results$msr, 
                  file.path(output_dir, "processed", "msr", 
                            paste0(base_name, "_msr.tiff")),
                  bits.per.sample = 16,
                  compression = "none",
                  native = TRUE)
      }
      
      # Bilateral
      if ("bilateral" %in% selected_methods) {
        cat("Applying bilateral filter...\n")
        results$bilateral <- matrix(apply_bilateral(curr_img),
                                    nrow = nrow(curr_img), 
                                    ncol = ncol(curr_img))
        writeTIFF(results$bilateral, 
                  file.path(output_dir, "processed", "bilateral", 
                            paste0(base_name, "_bilateral.tiff")),
                  bits.per.sample = 16,
                  compression = "none",
                  native = TRUE)
      }
      
      # LCEP
      if ("lcep" %in% selected_methods) {
        results$lcep <- matrix(apply_lcep(curr_img),
                               nrow = nrow(curr_img), 
                               ncol = ncol(curr_img))
        writeTIFF(results$lcep, 
                  file.path(output_dir, "processed", "lcep", 
                            paste0(base_name, "_lcep.tiff")),
                  bits.per.sample = 16,
                  compression = "none",
                  native = TRUE)
      }
      
      all_results[[i]] <- results
      
    }, error = function(e) {
      cat("Error processing", basename(img_path), ":", e$message, "\n")
      print(e)
    })
  }
  
  # Create median value plot
  png(file.path(output_dir, "comparisons", "median_values.png"),
      width = 800, height = 400)
  
  plot(median_values$frame, median_values$original, 
       type = "l", col = "blue", 
       xlab = "Frame Number", ylab = "Median Temperature",
       main = "Median Temperature Values Across Frames",
       ylim = range(c(median_values$original, median_values$matched)))
  
  lines(median_values$frame, median_values$matched, col = "red")
  
  legend("topright", 
         legend = c("Original", "Histogram Matched"),
         col = c("blue", "red"),
         lty = 1)
  
  dev.off()
  
  # Create visualization
  n_images <- length(tiff_files)
  n_methods <- length(selected_methods)
  
  png(file.path(output_dir, "comparisons", "methods_comparison.png"),
      width = 300 * n_methods, height = 400 * n_images)
  
  layout_matrix <- matrix(1:(n_methods * n_images), 
                          nrow = n_images, 
                          ncol = n_methods, 
                          byrow = TRUE)
  layout(layout_matrix)
  par(mar = c(2, 2, 3, 2))
  
  for (i in seq_along(all_results)) {
    results <- all_results[[i]]
    
    for (method in selected_methods) {
      plot(as.raster(colorize_thermal(results[[method]])),
           main = if(i == 1) method else "", axes = FALSE)
    }
  }
  
  dev.off()
  
  cat("\nProcessing complete. Selected methods:", paste(selected_methods, collapse = ", "), "\n")
  cat("Processed images saved in:", file.path(output_dir, "processed"), "\n")
  cat("Comparison visualization saved in:", file.path(output_dir, "comparisons"), "\n")
}

# Example usage:
# All methods:
# process_thermal_images("./images", "./processed_images")

# Selected methods:
# process_thermal_images("./images", "./processed_images", 
#                       methods = c("original", "hist_matched", "clahe"))

# Single method:
# process_thermal_images("./images", "./processed_images", 
#                       methods = "bilateral")