library(tiff)
library(EBImage)
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
  
  if (!dir.exists(input_dir)) stop("Input directory does not exist")
  subdirs <- c("processed", "comparisons")
  sapply(file.path(output_dir, subdirs), dir.create, recursive = TRUE, showWarnings = FALSE)
  
  tiff_files <- sort(list.files(input_dir, pattern = "\\.tiff?$", full.names = TRUE, ignore.case = TRUE))
  if (length(tiff_files) == 0) stop("No TIFF files found")
  
  # Calculate global min/max
  all_mins <- all_maxs <- numeric(length(tiff_files))
  for (i in seq_along(tiff_files)) {
    img <- readTIFF(tiff_files[i], all = TRUE)[[1]]
    all_mins[i] <- min(img, na.rm = TRUE)
    all_maxs[i] <- max(img, na.rm = TRUE)
  }
  global_min <- min(all_mins)
  global_max <- max(all_maxs)
  
  colorize_thermal <- function(img) {
    normalized <- (img - global_min) / (global_max - global_min)
    normalized[normalized < 0] <- 0
    normalized[normalized > 1] <- 1
    colors <- colorRampPalette(c("black", "purple", "red", "yellow", "white"))(256)
    matrix(colors[floor(normalized * 255) + 1], nrow = nrow(img))
  }
  
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
    
    # Create spatial Gaussian kernel
    spatial_kernel <- matrix(0, window_size, window_size)
    for(i in 1:window_size) {
      for(j in 1:window_size) {
        x <- i - half_size - 1
        y <- j - half_size - 1
        spatial_kernel[i,j] <- exp(-(x^2 + y^2)/(2 * sigma_s^2))
      }
    }
    spatial_kernel <- spatial_kernel / sum(spatial_kernel)
    
    # Apply bilateral filter
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
  
  # Process each image with selected methods
  all_results <- list()
  reference_img <- NULL
  hist_breaks <- seq(global_min, global_max, length.out = 100)
  
  for (i in seq_along(tiff_files)) {
    img_path <- tiff_files[i]
    cat("Processing:", basename(img_path), "\n")
    
    tryCatch({
      curr_img <- readTIFF(img_path, all = TRUE)[[1]]
      if (any(is.na(curr_img))) curr_img[is.na(curr_img)] <- global_min
      
      results <- list()
      results$original <- curr_img
      
      # Sequential Histogram Matching
      if ("hist_matched" %in% selected_methods) {
        if (is.null(reference_img)) {
          reference_img <- curr_img
          results$hist_matched <- curr_img
        } else {
          ref_hist <- hist(reference_img, breaks = hist_breaks, plot = FALSE)
          ref_cdf <- cumsum(ref_hist$counts) / sum(ref_hist$counts)
          
          src_hist <- hist(curr_img, breaks = hist_breaks, plot = FALSE)
          src_cdf <- cumsum(src_hist$counts) / sum(src_hist$counts)
          
          matched_img <- curr_img
          for (j in 1:length(curr_img)) {
            bin <- findInterval(curr_img[j], hist_breaks)
            if (bin > 0 && bin <= length(src_cdf)) {
              target_cdf <- src_cdf[bin]
              new_bin <- which.min(abs(ref_cdf - target_cdf))
              matched_img[j] <- hist_breaks[new_bin]
            }
          }
          results$hist_matched <- matched_img
          reference_img <- matched_img
        }
      }
      
      # CLAHE
      if ("clahe" %in% selected_methods) {
        scaled_img <- (curr_img - global_min) / (global_max - global_min)
        img_eb <- Image(scaled_img)
        clahe_img <- clahe(img_eb)
        results$clahe <- as.array(clahe_img) * (global_max - global_min) + global_min
      }
      
      # MSR
      if ("msr" %in% selected_methods) {
        results$msr <- apply_msr(curr_img)
      }
      
      # Bilateral
      if ("bilateral" %in% selected_methods) {
        cat("Applying bilateral filter...\n")
        results$bilateral <- apply_bilateral(curr_img)
      }
      
      # LCEP
      if ("lcep" %in% selected_methods) {
        results$lcep <- apply_lcep(curr_img)
      }
      
      all_results[[i]] <- results
      
    }, error = function(e) {
      cat("Error processing", basename(img_path), ":", e$message, "\n")
      print(e)
    })
  }
  
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
}


process_thermal_images("./images", "./processed_images")

# Example usage:
# All methods:
# process_thermal_images("./images", "./processed_images")

# Selected methods:
# process_thermal_images("./images", "./processed_images", 
#                       methods = c("original", "hist_matched", "clahe"))

# Single method:
# process_thermal_images("./images", "./processed_images", 
#                       methods = "bilateral")