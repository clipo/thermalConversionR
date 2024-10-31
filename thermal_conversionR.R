# Required libraries
library(tiff)
library(EBImage)
library(grid)
library(gridExtra)
library(parallel)
library(doParallel)
library(stats)

process_thermal_images <- function(input_dir, output_dir, methods = "all", parallel = TRUE) {
  # Start time tracking
  start_time <- Sys.time()
  
  # Available methods - Expanded list
  available_methods <- c(
    "original", 
    "hist_matched", 
    "clahe", 
    "msr", 
    "bilateral", 
    "lcep",
    "adaptive_gamma",
    "guided_filter",
    "wavelet_fusion",
    "local_histogram"
  )
  
  # Setup parallel processing if requested
  if (parallel) {
    num_cores <- parallel::detectCores() - 1
    if (num_cores > 1) {
      cl <- makeCluster(num_cores)
      registerDoParallel(cl)
      cat("Parallel processing enabled with", num_cores, "cores\n")
    }
  }
  
  # Process methods argument
  if (length(methods) == 1 && methods == "all") {
    selected_methods <- available_methods
  } else {
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
  
  # Directory validation and setup
  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }
  
  # Get list of TIFF files with error handling
  tiff_files <- list.files(input_dir, 
                           pattern = "\\.tiff?$", 
                           full.names = TRUE, 
                           ignore.case = TRUE)
  
  if (length(tiff_files) == 0) {
    stop("No TIFF files found in directory: ", input_dir)
  }
  
  # Create output directory structure
  subdirs <- c("comparisons", "processed", "metrics")
  for (dir in subdirs) {
    dir.create(file.path(output_dir, dir), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create method-specific subdirectories
  method_dirs <- setdiff(selected_methods, "original")
  sapply(file.path(output_dir, "processed", method_dirs), 
         dir.create, recursive = TRUE, showWarnings = FALSE)
  
  # Helper functions
  read_tiff_safe <- function(path) {
    tryCatch({
      img <- readTIFF(path, all = TRUE)[[1]]
      gc()  # Force garbage collection
      return(img)
    }, error = function(e) {
      stop(paste("Error reading file:", path, "-", e$message))
    })
  }
  
  match_histogram <- function(source_img, reference_img) {
    # Calculate histograms
    hist_breaks <- seq(global_min, global_max, length.out = 256)
    
    # Reference histogram
    ref_hist <- hist(reference_img, breaks = hist_breaks, plot = FALSE)
    ref_cdf <- cumsum(ref_hist$counts) / sum(ref_hist$counts)
    
    # Source histogram
    src_hist <- hist(source_img, breaks = hist_breaks, plot = FALSE)
    src_cdf <- cumsum(src_hist$counts) / sum(src_hist$counts)
    
    # Create lookup table for matching
    lookup <- numeric(length(hist_breaks) - 1)
    for(i in seq_along(lookup)) {
      # Find the closest CDF value in reference
      lookup[i] <- hist_breaks[which.min(abs(ref_cdf - src_cdf[i]))]
    }
    
    # Apply lookup table to source image
    result <- matrix(source_img, nrow = nrow(source_img), ncol = ncol(source_img))
    for(i in seq_along(source_img)) {
      bin <- findInterval(source_img[i], hist_breaks)
      if(bin > 0 && bin <= length(lookup)) {
        result[i] <- lookup[bin]
      }
    }
    
    return(result)
  }
  
  colorize_thermal <- function(img, colormap = "inferno") {
    normalized <- (img - global_min) / (global_max - global_min)
    normalized[normalized < 0] <- 0
    normalized[normalized > 1] <- 1
    
    colormaps <- list(
      inferno = c("black", "purple", "red", "yellow", "white"),
      viridis = c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725"),
      magma = c("black", "#221150", "#B63679", "#FB8861", "#FCFDBF"),
      plasma = c("#0D0887", "#6A00A8", "#B12A90", "#E16462", "#FCA636", "#F0F921")
    )
    
    if (!(colormap %in% names(colormaps))) {
      warning("Invalid colormap, using inferno")
      colormap <- "inferno"
    }
    
    colors <- colorRampPalette(colormaps[[colormap]])(256)
    matrix(colors[floor(normalized * 255) + 1], nrow = nrow(img))
  }
  
  # Enhanced MSR function
  apply_msr <- function(img, scales = c(2, 4, 8), weights = NULL) {
    if (is.null(weights)) {
      weights <- rep(1/length(scales), length(scales))
    }
    
    result <- matrix(0, nrow = nrow(img), ncol = ncol(img))
    img_norm <- (img - min(img)) / (max(img) - min(img))
    img_eb <- Image(img_norm)
    
    for (i in seq_along(scales)) {
      gaussian <- gblur(img_eb, sigma = scales[i])
      result <- result + weights[i] * (log(img_norm + 1) - log(as.matrix(gaussian) + 1))
    }
    
    result <- (result - min(result)) / (max(result) - min(result))
    lower <- quantile(result, 0.01)
    upper <- quantile(result, 0.99)
    result <- (result - lower) / (upper - lower)
    result[result < 0] <- 0
    result[result > 1] <- 1
    
    result <- result * (global_max - global_min) + global_min
    return(result)
  }
  
  # Enhanced bilateral filter
  apply_bilateral <- function(img, sigma_s = NULL, sigma_r = NULL) {
    if (is.null(sigma_s)) {
      sigma_s <- min(nrow(img), ncol(img)) / 50
    }
    if (is.null(sigma_r)) {
      sigma_r <- (max(img) - min(img)) / 10
    }
    
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
    
    cat("\nApplying bilateral filter...\n")
    pb <- txtProgressBar(min = 0, max = nrow(img_norm), style = 3)
    
    for(i in 1:nrow(img_norm)) {
      setTxtProgressBar(pb, i)
      for(j in 1:ncol(img_norm)) {
        r1 <- max(1, i - half_size)
        r2 <- min(nrow(img_norm), i + half_size)
        c1 <- max(1, j - half_size)
        c2 <- min(ncol(img_norm), j + half_size)
        
        window <- img_norm[r1:r2, c1:c2]
        intensity_diff <- (window - img_norm[i,j])^2
        range_kernel <- exp(-intensity_diff/(2 * sigma_r^2))
        weights <- spatial_kernel[1:(r2-r1+1), 1:(c2-c1+1)] * range_kernel
        result[i,j] <- sum(window * weights) / sum(weights)
      }
    }
    close(pb)
    
    result <- result * (max(img) - min(img)) + min(img)
    return(result)
  }
  
  # Additional processing methods
  apply_adaptive_gamma <- function(img, window_size = 51) {
    # Normalize input image
    img_norm <- (img - min(img)) / (max(img) - min(img))
    
    cat("\nApplying adaptive gamma correction...\n")
    
    tryCatch({
      # Create padded image for convolution
      pad_size <- (window_size - 1) %/% 2
      padded_img <- matrix(0, 
                           nrow = nrow(img_norm) + 2 * pad_size,
                           ncol = ncol(img_norm) + 2 * pad_size)
      padded_img[(pad_size + 1):(pad_size + nrow(img_norm)),
                 (pad_size + 1):(pad_size + ncol(img_norm))] <- img_norm
      
      # Calculate local statistics manually
      cat("Calculating local statistics...\n")
      local_mean <- matrix(0, nrow = nrow(img_norm), ncol = ncol(img_norm))
      local_std <- matrix(0, nrow = nrow(img_norm), ncol = ncol(img_norm))
      
      # Set up progress bar
      pb <- txtProgressBar(min = 0, max = nrow(img_norm), style = 3)
      
      # Calculate local statistics using sliding window
      for(i in 1:nrow(img_norm)) {
        for(j in 1:ncol(img_norm)) {
          # Extract window
          window <- padded_img[i:(i + window_size - 1),
                               j:(j + window_size - 1)]
          
          # Calculate statistics
          local_mean[i,j] <- mean(window)
          local_std[i,j] <- sd(window)
        }
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      cat("\nComputing adaptive gamma map...\n")
      
      # Compute adaptive gamma values
      gamma_map <- 1 + (0.5 * (1 - local_mean)) * (1 + local_std)
      
      # Limit gamma range
      gamma_map <- pmin(pmax(gamma_map, 0.5), 2.0)
      
      cat("Applying gamma correction...\n")
      
      # Apply gamma correction
      result <- img_norm^gamma_map
      
      # Scale back to original range
      result <- result * (global_max - global_min) + global_min
      
      cat("Adaptive gamma correction completed\n")
      return(result)
      
    }, error = function(e) {
      cat("\nError in adaptive gamma correction:\n")
      cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
      cat("Window size:", window_size, "\n")
      cat("Error message:", e$message, "\n")
      stop(e)
    })
  }
  
  apply_guided_filter <- function(img, radius = NULL, epsilon = 1e-6) {
    # Ensure radius is valid and odd
    if (is.null(radius)) {
      radius <- min(nrow(img), ncol(img)) %/% 50
    }
    # Make sure radius is odd
    radius <- if(radius %% 2 == 0) radius + 1 else radius
    
    cat(sprintf("\nUsing radius: %d\n", radius))
    
    img_norm <- (img - min(img)) / (max(img) - min(img))
    guide <- img_norm
    
    cat("Applying guided filter...\n")
    
    tryCatch({
      # Manual box filter implementation
      box_filter <- function(input, r) {
        # Add padding
        h <- nrow(input)
        w <- ncol(input)
        pad_size <- r %/% 2
        
        padded <- matrix(0, h + 2*pad_size, w + 2*pad_size)
        padded[(pad_size+1):(pad_size+h), (pad_size+1):(pad_size+w)] <- input
        
        output <- matrix(0, h, w)
        window_area <- (2*pad_size + 1)^2
        
        cat("Filtering... \n")
        pb <- txtProgressBar(min = 0, max = h, style = 3)
        
        for(i in 1:h) {
          for(j in 1:w) {
            window <- padded[i:(i+2*pad_size), j:(j+2*pad_size)]
            output[i,j] <- sum(window) / window_area
          }
          setTxtProgressBar(pb, i)
        }
        close(pb)
        
        return(output)
      }
      
      # Calculate means using box filter
      mean_I <- box_filter(guide, radius)
      mean_p <- box_filter(img_norm, radius)
      mean_II <- box_filter(guide * guide, radius)
      mean_Ip <- box_filter(guide * img_norm, radius)
      
      # Calculate variance and covariance
      var_I <- mean_II - mean_I * mean_I
      cov_Ip <- mean_Ip - mean_I * mean_p
      
      # Calculate coefficients
      a <- cov_Ip / (var_I + epsilon)
      b <- mean_p - a * mean_I
      
      # Filter coefficients
      mean_a <- box_filter(a, radius)
      mean_b <- box_filter(b, radius)
      
      # Final guided filter result
      result <- mean_a * guide + mean_b
      
      # Scale back to original range
      result <- result * (global_max - global_min) + global_min
      
      cat("\nGuided filter completed\n")
      return(result)
      
    }, error = function(e) {
      cat("\nError in guided filter:\n")
      cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
      cat("Radius:", radius, "\n")
      cat("Error message:", e$message, "\n")
      stop(e)
    })
  }
  
  apply_wavelet_fusion <- function(img, levels = 3, enhancement_factor = 1.2) {
    img_norm <- (img - min(img)) / (max(img) - min(img))
    
    cat("\nPerforming multi-scale decomposition...\n")
    
    tryCatch({
      # Simplified Gaussian blur with better boundary handling
      gaussian_blur <- function(img, sigma = 1.5) {
        kernel_size <- max(3, ceiling(6 * sigma))
        if(kernel_size %% 2 == 0) kernel_size <- kernel_size + 1
        
        x <- seq(-(kernel_size%/%2), kernel_size%/%2)
        kernel <- dnorm(x, 0, sigma)
        kernel <- kernel / sum(kernel)
        
        h <- nrow(img)
        w <- ncol(img)
        result <- img
        
        # Horizontal blur
        for(i in 1:h) {
          row <- img[i,]
          # Replicate border for padding
          padded_row <- c(rep(row[1], kernel_size%/%2), 
                          row, 
                          rep(row[w], kernel_size%/%2))
          
          # Convolve
          for(j in 1:w) {
            result[i,j] <- sum(padded_row[j:(j+kernel_size-1)] * kernel)
          }
        }
        
        # Vertical blur
        temp <- result
        for(j in 1:w) {
          col <- temp[,j]
          # Replicate border for padding
          padded_col <- c(rep(col[1], kernel_size%/%2), 
                          col, 
                          rep(col[h], kernel_size%/%2))
          
          # Convolve
          for(i in 1:h) {
            result[i,j] <- sum(padded_col[i:(i+kernel_size-1)] * kernel)
          }
        }
        
        return(result)
      }
      
      # Create Gaussian pyramid with dimension checks
      pyramids <- list()
      current <- img_norm
      pyramids[[1]] <- current
      
      cat("Building pyramid levels...\n")
      pb <- txtProgressBar(min = 0, max = levels, style = 3)
      
      for(i in 2:levels) {
        # Check if image is too small
        if(min(dim(current)) < 4) {
          warning("Image too small for further decomposition. Stopping at level ", i-1)
          levels <- i-1
          break
        }
        
        # Gaussian blur and downsample
        blurred <- gaussian_blur(current)
        current <- blurred[seq(1, nrow(blurred), 2), seq(1, ncol(blurred), 2)]
        pyramids[[i]] <- current
        
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      # Process each level
      cat("\nEnhancing pyramid levels...\n")
      enhanced_pyramids <- list()
      pb <- txtProgressBar(min = 0, max = levels, style = 3)
      
      for(i in 1:levels) {
        current_level <- pyramids[[i]]
        local_mean <- gaussian_blur(current_level)
        detail <- current_level - local_mean
        
        # Adaptive enhancement
        level_enhancement <- enhancement_factor * (1 + 0.5 * (levels - i) / levels)
        enhanced_detail <- detail * level_enhancement
        
        enhanced_pyramids[[i]] <- local_mean + enhanced_detail
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      # Reconstruct image
      cat("\nReconstructing final image...\n")
      result <- enhanced_pyramids[[levels]]
      
      for(i in (levels-1):1) {
        # Simple upsampling with boundary checking
        current_small <- result
        small_h <- nrow(current_small)
        small_w <- ncol(current_small)
        
        # Target dimensions from the current pyramid level
        target_h <- nrow(enhanced_pyramids[[i]])
        target_w <- ncol(enhanced_pyramids[[i]])
        
        # Nearest neighbor upsampling (more stable than bilinear)
        upsampled <- matrix(0, target_h, target_w)
        for(y in 1:target_h) {
          src_y <- min(max(1, ceiling(y/2)), small_h)
          for(x in 1:target_w) {
            src_x <- min(max(1, ceiling(x/2)), small_w)
            upsampled[y,x] <- current_small[src_y, src_x]
          }
        }
        
        # Smooth the upsampled result
        upsampled_smooth <- gaussian_blur(upsampled, sigma = 0.8)
        
        # Add the details from the current level
        result <- enhanced_pyramids[[i]] + upsampled_smooth
      }
      
      # Final normalization with boundary checking
      result[result < 0] <- 0
      result[result > 1] <- 1
      result <- result * (global_max - global_min) + global_min
      
      cat("Multi-scale enhancement completed\n")
      return(result)
      
    }, error = function(e) {
      cat("\nError in multi-scale enhancement:\n")
      cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
      cat("Number of levels:", levels, "\n")
      cat("Error message:", e$message, "\n")
      cat("Stack trace:\n")
      print(sys.calls())
      stop(e)
    })
  }
  
  # Add this helper function inside process_thermal_images
  apply_clahe <- function(img, num_tiles = NULL, clip_limit = 0.03, smooth_factor = 3) {
    # Memory efficient Gaussian blur
    gaussian_blur <- function(img, sigma) {
      # Separate 1D kernels for memory efficiency
      kernel_size <- max(3, ceiling(sigma * 6))
      if(kernel_size %% 2 == 0) kernel_size <- kernel_size + 1
      
      x <- seq(-kernel_size/2, kernel_size/2)
      kernel_1d <- dnorm(x, 0, sigma)
      kernel_1d <- kernel_1d / sum(kernel_1d)
      
      # Horizontal pass
      temp <- matrix(0, nrow(img), ncol(img))
      pad <- kernel_size %/% 2
      
      # Process rows
      for(i in 1:nrow(img)) {
        # Pad row by replication
        row <- img[i,]
        padded_row <- c(rep(row[1], pad), row, rep(row[length(row)], pad))
        
        # Convolve
        for(j in 1:ncol(img)) {
          temp[i,j] <- sum(padded_row[j:(j+kernel_size-1)] * kernel_1d)
        }
      }
      
      # Vertical pass
      result <- matrix(0, nrow(img), ncol(img))
      
      # Process columns
      for(j in 1:ncol(temp)) {
        # Pad column by replication
        col <- temp[,j]
        padded_col <- c(rep(col[1], pad), col, rep(col[length(col)], pad))
        
        # Convolve
        for(i in 1:nrow(img)) {
          result[i,j] <- sum(padded_col[i:(i+kernel_size-1)] * kernel_1d)
        }
      }
      
      return(result)
    }
    
    # Calculate adaptive number of tiles
    if (is.null(num_tiles)) {
      suggested_tiles <- min(max(floor(min(dim(img)) / 128), 3), 8)  # Reduced max tiles
      num_tiles_y <- suggested_tiles
      num_tiles_x <- suggested_tiles
    } else {
      num_tiles_y <- num_tiles
      num_tiles_x <- num_tiles
    }
    
    img_norm <- (img - min(img)) / (max(img) - min(img))
    height <- nrow(img)
    width <- ncol(img)
    
    # Process tiles with overlap
    result <- matrix(0, height, width)
    weights <- matrix(0, height, width)
    
    cat(sprintf("\nApplying CLAHE (tiles: %dx%d, clip limit: %.3f)...\n", 
                num_tiles_x, num_tiles_y, clip_limit))
    
    pb <- txtProgressBar(min = 0, max = num_tiles_y, style = 3)
    
    # Process each tile
    for(ty in 1:num_tiles_y) {
      setTxtProgressBar(pb, ty)
      
      # Calculate tile boundaries with overlap
      y_start <- max(1, floor((ty-1) * height/num_tiles_y - height/(num_tiles_y * 4)))
      y_end <- min(height, ceiling(ty * height/num_tiles_y + height/(num_tiles_y * 4)))
      
      for(tx in 1:num_tiles_x) {
        x_start <- max(1, floor((tx-1) * width/num_tiles_x - width/(num_tiles_x * 4)))
        x_end <- min(width, ceiling(tx * width/num_tiles_x + width/(num_tiles_x * 4)))
        
        # Extract and process tile
        tile <- img_norm[y_start:y_end, x_start:x_end]
        
        # Calculate histogram
        hist_result <- hist(tile, breaks = seq(0, 1, length.out = 257), plot = FALSE)
        hist_counts <- hist_result$counts
        
        # Apply clipping
        clip_height <- clip_limit * prod(dim(tile))
        excess <- sum(pmax(hist_counts - clip_height, 0))
        redistrib_amt <- excess / length(hist_counts)
        
        hist_limited <- pmin(hist_counts, clip_height) + redistrib_amt
        
        # Calculate CDF
        cdf <- cumsum(hist_limited)
        cdf <- (cdf - min(cdf)) / (max(cdf) - min(cdf))
        
        # Create weight matrix (Gaussian weights for smooth blending)
        y_rel <- seq(0, 1, length.out = y_end - y_start + 1)
        x_rel <- seq(0, 1, length.out = x_end - x_start + 1)
        y_weights <- dnorm(y_rel, mean = 0.5, sd = 0.25)
        x_weights <- dnorm(x_rel, mean = 0.5, sd = 0.25)
        weight_matrix <- outer(y_weights, x_weights)
        weight_matrix <- weight_matrix / max(weight_matrix)
        
        # Transform tile
        tile_result <- matrix(0, nrow(tile), ncol(tile))
        for(i in 1:nrow(tile)) {
          val <- tile[i,]
          bins <- pmin(ceiling(val * 255) + 1, 256)
          tile_result[i,] <- cdf[bins]
        }
        
        # Accumulate weighted results
        result[y_start:y_end, x_start:x_end] <- 
          result[y_start:y_end, x_start:x_end] + tile_result * weight_matrix
        weights[y_start:y_end, x_start:x_end] <- 
          weights[y_start:y_end, x_start:x_end] + weight_matrix
      }
    }
    close(pb)
    
    # Normalize by weights
    result <- result / (weights + .Machine$double.eps)
    
    # Apply smoothing passes
    cat("\nApplying smoothing passes...\n")
    smooth_result <- result
    if (smooth_factor>1) {
      for(i in 1:smooth_factor) {
        sigma <- min(width, height) / (50 * sqrt(i))  # Decreasing sigma for each pass
        smooth_result <- gaussian_blur(smooth_result, sigma)
        cat(sprintf("Smoothing pass %d/%d completed\n", i, smooth_factor))
        gc()  # Force garbage collection
      }
    }
    # Scale back to original range
    smooth_result <- smooth_result * (global_max - global_min) + global_min
    
    cat("CLAHE completed\n")
    return(smooth_result)
  }
  
  apply_local_histogram <- function(img, window_size = NULL) {
    # Determine window size if not provided
    if (is.null(window_size)) {
      window_size <- min(nrow(img), ncol(img)) %/% 20
      if (window_size %% 2 == 0) window_size <- window_size + 1
      window_size <- max(3, window_size)  # Ensure minimum size of 3
    }
    
    cat(sprintf("\nApplying local histogram equalization (window size: %d)...\n", window_size))
    
    tryCatch({
      # Normalize input image
      img_norm <- (img - min(img)) / (max(img) - min(img))
      result <- matrix(0, nrow = nrow(img), ncol = ncol(img))
      
      # Calculate padding size
      pad_size <- window_size %/% 2
      
      # Create padded image using border replication
      padded <- matrix(0, 
                       nrow = nrow(img) + 2 * pad_size,
                       ncol = ncol(img) + 2 * pad_size)
      
      # Fill center with actual image
      padded[(pad_size + 1):(pad_size + nrow(img)),
             (pad_size + 1):(pad_size + ncol(img))] <- img_norm
      
      # Replicate borders
      # Top and bottom padding
      padded[1:pad_size, (pad_size + 1):(pad_size + ncol(img))] <- 
        matrix(img_norm[1,], pad_size, ncol(img), byrow = TRUE)
      padded[(nrow(padded) - pad_size + 1):nrow(padded), 
             (pad_size + 1):(pad_size + ncol(img))] <- 
        matrix(img_norm[nrow(img),], pad_size, ncol(img), byrow = TRUE)
      
      # Left and right padding without pipes
      left_border <- matrix(padded[, pad_size + 1], nrow = nrow(padded), ncol = pad_size)
      right_border <- matrix(padded[, ncol(padded) - pad_size], 
                             nrow = nrow(padded), ncol = pad_size)
      padded[, 1:pad_size] <- left_border
      padded[, (ncol(padded) - pad_size + 1):ncol(padded)] <- right_border
      
      # Process image with progress bar
      cat("Processing image...\n")
      pb <- txtProgressBar(min = 0, max = nrow(img), style = 3)
      
      # Local histogram equalization
      for(i in 1:nrow(img)) {
        for(j in 1:ncol(img)) {
          # Extract local window
          window <- padded[i:(i + window_size - 1),
                           j:(j + window_size - 1)]
          
          # Calculate local CDF
          hist_counts <- hist(window, breaks = 256, plot = FALSE)$counts
          cdf <- cumsum(hist_counts) / sum(hist_counts)
          
          # Transform pixel using local CDF
          pixel_bin <- ceiling(img_norm[i,j] * 255) + 1
          result[i,j] <- cdf[min(pixel_bin, 256)]
        }
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      # Apply contrast limiting to reduce noise amplification
      clip_limit <- 0.03  # Adjust this value to control contrast limitation
      result_hist <- hist(result, breaks = 256, plot = FALSE)$counts
      excess <- pmax(result_hist - (clip_limit * length(result)), 0)
      redistribution <- sum(excess) / 256
      
      # Create lookup table for contrast limiting
      hist_clipped <- pmin(result_hist, clip_limit * length(result)) + redistribution
      cdf_clipped <- cumsum(hist_clipped) / sum(hist_clipped)
      
      # Apply contrast limiting
      for(i in 1:length(result)) {
        bin <- ceiling(result[i] * 255) + 1
        result[i] <- cdf_clipped[min(bin, 256)]
      }
      
      # Scale back to original range
      result <- result * (global_max - global_min) + global_min
      
      cat("\nLocal histogram equalization completed\n")
      return(result)
      
    }, error = function(e) {
      cat("\nError in local histogram equalization:\n")
      cat("Image dimensions:", nrow(img), "x", ncol(img), "\n")
      cat("Window size:", window_size, "\n")
      cat("Error message:", e$message, "\n")
      stop(e)
    })
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
  
  create_comparison_visualization <- function(all_results, selected_methods, tiff_files, output_dir, max_images = 100) {
    # Limit the number of images to process
    n_images <- min(length(all_results), max_images)
    if (length(all_results) > max_images) {
      cat(sprintf("\nWarning: Limiting visualization to first %d images out of %d total images\n", 
                  max_images, length(all_results)))
    }
    
    n_methods <- length(selected_methods)
    
    # Calculate dimensions with minimum sizes
    plot_width <- max(800, 300 * n_methods)  # Minimum width of 800 pixels
    plot_height <- max(600, 200 * n_images + 100)  # Minimum height of 600 pixels
    
    cat(sprintf("\nCreating comparison visualization (%d images x %d methods)...\n", 
                n_images, n_methods))
    
    # Create the PNG device with adjusted dimensions
    png_file <- file.path(output_dir, "comparisons", "methods_comparison.png")
    png(png_file,
        width = plot_width, 
        height = plot_height,
        res = 100)  # Reduced resolution for better margin handling
    
    # Create layout for images
    layout_matrix <- matrix(1:(n_methods * n_images), 
                            nrow = n_images, 
                            ncol = n_methods, 
                            byrow = TRUE)
    layout(layout_matrix)
    
    # Create color scale for temperature mapping
    temp_colors <- colorRampPalette(c("black", "purple", "red", "yellow", "white"))(100)
    
    # Process only the first max_images images
    for(i in 1:n_images) {
      results <- all_results[[i]]
      file_name <- basename(tiff_files[i])
      
      for(method in selected_methods) {
        # Adjust margins based on position
        if(method == selected_methods[1]) {
          # Left column - space for filename
          par(mar = c(1, 2, 3, 1))
        } else if(method == tail(selected_methods, 1)) {
          # Right column - space for temperature scale
          par(mar = c(1, 1, 3, 3))
        } else {
          # Middle columns - minimal margins
          par(mar = c(1, 1, 3, 1))
        }
        
        # Plot image
        img_data <- results[[method]]
        
        # Create title
        if(i == 1) {
          main_title <- method
        } else {
          main_title <- ""
        }
        
        if(method == selected_methods[1]) {
          main_title <- paste0(main_title, "\n", 
                               sprintf("(%d/%d) %s", i, n_images, file_name))
        }
        
        # Plot with enhanced temperature scale
        image(1:ncol(img_data), 1:nrow(img_data), 
              t(img_data[nrow(img_data):1,]), 
              col = temp_colors,
              main = main_title,
              xlab = "", ylab = "",
              axes = FALSE)
        
        # Add temperature scale bar for rightmost column
        if(method == tail(selected_methods, 1)) {
          temp_range <- pretty(c(global_min, global_max), n = 5)
          axis(4, at = seq(0, 1, length.out = length(temp_range)),
               labels = sprintf("%.1f°C", temp_range),
               las = 2,
               cex.axis = 0.8)  # Smaller axis text
        }
      }
      
      # Update progress
      if (i %% 10 == 0) {
        cat(sprintf("Processed %d/%d images\n", i, n_images))
      }
    }
    
    # Close the device
    dev.off()
    
    cat(sprintf("\nVisualization complete. Saved to: %s\n", png_file))
    
    # Create a separate summary file
    summary_file <- file.path(output_dir, "comparisons", "processing_summary.txt")
    cat(sprintf("Creating summary file: %s\n", summary_file))
    
    writeLines(
      c(sprintf("Thermal Image Processing Summary"),
        sprintf("============================"),
        sprintf("Total images processed: %d", length(all_results)),
        sprintf("Images in visualization: %d", n_images),
        sprintf("Methods applied: %s", paste(selected_methods, collapse = ", ")),
        sprintf("\nProcessed Files:"),
        sprintf("---------------"),
        sapply(seq_along(tiff_files), function(i) {
          sprintf("%3d. %s %s", 
                  i, 
                  basename(tiff_files[i]),
                  if(i <= max_images) "(included in visualization)" else "(not visualized)")
        })),
      summary_file
    )
  }
  
  # Calculate global statistics
  cat("\nCalculating global statistics...\n")
  pb <- txtProgressBar(min = 0, max = length(tiff_files), style = 3)
  
  global_stats <- lapply(seq_along(tiff_files), function(i) {
    img <- read_tiff_safe(tiff_files[i])
    setTxtProgressBar(pb, i)
    c(min = min(img, na.rm = TRUE), 
      max = max(img, na.rm = TRUE),
      mean = mean(img, na.rm = TRUE),
      sd = sd(img, na.rm = TRUE))
  })
  close(pb)
  
  global_min <- min(sapply(global_stats, function(x) x["min"]))
  global_max <- max(sapply(global_stats, function(x) x["max"]))
  global_mean <- mean(sapply(global_stats, function(x) x["mean"]))
  global_sd <- mean(sapply(global_stats, function(x) x["sd"]))
  
  # Main processing function
  process_images <- function() {
    all_results <- list()
    
    for(i in seq_along(tiff_files)) {
      img_path <- tiff_files[i]
      base_name <- tools::file_path_sans_ext(basename(img_path))
      cat("\nProcessing image", i, "of", length(tiff_files), ":", base_name, "\n")
      
      tryCatch({
        curr_img <- read_tiff_safe(img_path)
        if (any(is.na(curr_img))) curr_img[is.na(curr_img)] <- global_min
        
        results <- list()
        results$original <- curr_img
        
        for(method in setdiff(selected_methods, "original")) {
          cat("Applying", method, "...\n")
          results[[method]] <- switch(method,
                                      "hist_matched" = {
                                        if(i == 1) {
                                          reference_img <<- curr_img
                                          curr_img
                                        } else {
                                          match_histogram(curr_img, reference_img)
                                        }
                                      },
                                      "clahe" = {
                                        apply_clahe(curr_img,
                                                    num_tiles = NULL,     # Automatic (3-8 tiles)
                                                    clip_limit = 0.03,    # Contrast limit
                                                    smooth_factor = 0)    # Number of smoothing passes
                                      },
                                      "msr" = apply_msr(curr_img),
                                      "bilateral" = apply_bilateral(curr_img),
                                      "lcep" = apply_lcep(curr_img),
                                      "adaptive_gamma" = apply_adaptive_gamma(curr_img),
                                      "guided_filter" = apply_guided_filter(curr_img),
                                      "wavelet_fusion" = apply_wavelet_fusion(curr_img),
                                      "local_histogram" = apply_local_histogram(curr_img)
          )
          
          writeTIFF(results[[method]], 
                    file.path(output_dir, "processed", method, 
                              paste0(base_name, "_", method, ".tiff")),
                    bits.per.sample = 16,
                    compression = "none")
        }
        
        all_results[[i]] <- results
        
      }, error = function(e) {
        cat("Error processing", basename(img_path), ":", e$message, "\n")
        print(e)
      })
    }
    
    return(all_results)
  }
  
  # Run main processing
  all_results <- process_images()
  
  # Create visualizations
  create_comparison_visualization(all_results, selected_methods, tiff_files, output_dir)
  
  # Clean up parallel processing
  if (exists("cl")) {
    stopCluster(cl)
  }
  
  # Report processing time
  end_time <- Sys.time()
  processing_time <- difftime(end_time, start_time, units = "mins")
  cat("\nProcessing completed in", round(processing_time, 2), "minutes\n")
  
  return(invisible(all_results))
}