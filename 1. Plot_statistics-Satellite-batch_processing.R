# Install required packages if needed
# install.packages(c("terra", "sf", "dplyr", "ggplot2", "tidyr", "magrittr"))
library("terra")
library("sf")
library("dplyr")
library("ggplot2")
library("tidyr")
library("magrittr")

# Increase memory limits for better stability
options(terra.memfrac=0.5)  # Use 50% of available memory for terra operations
if (.Platform$OS.type == "windows") {
  memory.limit(size=8000)  # Adjust based on your system's RAM
}

#----------------------------------------------------------------
# PRINT DIRECTORY FUNCTION
#----------------------------------------------------------------
print_dir_three_levels <- function(base_dir = getwd()) {
  cat(basename(base_dir), "\n")
  level1_items <- list.files(base_dir, full.names = FALSE)
  
  for (item1 in sort(level1_items)) {
    item1_path <- file.path(base_dir, item1)
    cat("|-- ", item1, "\n", sep = "")
    
    if (file.info(item1_path)$isdir) {
      level2_items <- list.files(item1_path, full.names = FALSE)
      
      for (item2 in sort(level2_items)) {
        item2_path <- file.path(item1_path, item2)
        cat("|   |-- ", item2, "\n", sep = "")
        
        if (file.info(item2_path)$isdir) {
          level3_items <- list.files(item2_path, full.names = FALSE)
          
          for (item3 in sort(level3_items)) {
            cat("|   |   |-- ", item3, "\n", sep = "")
          }
        }
      }
    }
  }
}

#----------------------------------------------------------------
# PROCESS SINGLE SATELLITE FOLDER FUNCTION
#----------------------------------------------------------------
process_single_satellite_folder <- function(base_dir, field, date_folder_name) {
  # Create log file
  log_file <- file.path(base_dir, paste0("satellite_processing_log_", field, "_", date_folder_name, ".txt"))
  cat(paste("Processing started at:", Sys.time(), "\n"), file=log_file)
  
  # Construct path to specified folder
  field_path <- file.path(base_dir, field)
  if (!dir.exists(field_path)) {
    error_msg <- paste("Field directory not found:", field_path)
    message(error_msg)
    cat(paste(error_msg, "\n"), file=log_file, append=TRUE)
    return(FALSE)
  }
  
  # Path to Satellite folder
  satellite_path <- file.path(field_path, "Satellite")
  if (!dir.exists(satellite_path)) {
    error_msg <- paste("Satellite directory not found for field:", field)
    message(error_msg)
    cat(paste(error_msg, "\n"), file=log_file, append=TRUE)
    return(FALSE)
  }
  
  # Path to specific date folder
  date_folder <- file.path(satellite_path, date_folder_name)
  if (!dir.exists(date_folder)) {
    error_msg <- paste("Date folder not found:", date_folder)
    message(error_msg)
    cat(paste(error_msg, "\n"), file=log_file, append=TRUE)
    return(FALSE)
  }
  
  # Process the folder
  message(paste("Processing satellite folder:", date_folder_name))
  cat(paste("\nProcessing satellite folder:", date_folder_name, "at", Sys.time(), "\n"), file=log_file, append=TRUE)
  
  # Find the necessary files in this folder
  tryCatch({
    # Match pattern for satellite imagery (differs based on field)
    if (field == "SF") {
      satellite_files <- list.files(date_folder, pattern="AOI_Soyface.*\\.tif$", 
                                    full.names=TRUE, ignore.case=TRUE)
    } else if (field == "EF") {
      satellite_files <- list.files(date_folder, pattern="AOI_EF.*\\.tif$", 
                                    full.names=TRUE, ignore.case=TRUE)
    } else {
      satellite_files <- list.files(date_folder, pattern="AOI_.*\\.tif$", 
                                    full.names=TRUE, ignore.case=TRUE)
    }
    
    # Check for Grid_files subfolder for shapefiles
    grid_folder <- file.path(date_folder, "Grid_files")
    if (dir.exists(grid_folder)) {
      shp_files <- list.files(grid_folder, pattern="\\.shp$", full.names=TRUE)
    } else {
      # Try in the main folder
      shp_files <- list.files(date_folder, pattern="\\.shp$", full.names=TRUE)
    }
    
    # Check if we found all required files
    if (length(satellite_files) == 0) {
      error_msg <- paste("No satellite image file found in:", date_folder)
      message(error_msg)
      cat(paste(error_msg, "\n"), file=log_file, append=TRUE)
      return(FALSE)
    }
    if (length(shp_files) == 0) {
      error_msg <- paste("No shapefile found in:", date_folder, "or its Grid_files subfolder")
      message(error_msg)
      cat(paste(error_msg, "\n"), file=log_file, append=TRUE)
      return(FALSE)
    }
    
    # Use the first file of each type if multiple found
    satellite_path <- satellite_files[1]
    shapefile_path <- shp_files[1]
    
    # Create output folder inside the date folder
    output_folder <- file.path(date_folder, "results")
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive=TRUE)
    }
    
    # Log the files found
    cat(paste("Satellite image file:", basename(satellite_path), "\n"), file=log_file, append=TRUE)
    cat(paste("Shapefile:", basename(shapefile_path), "\n"), file=log_file, append=TRUE)
    cat(paste("Output folder:", output_folder, "\n"), file=log_file, append=TRUE)
    
    # Determine band configuration based on date
    date_str <- as.character(date_folder_name)
    # Extract the date part with format
    date_parts <- strsplit(date_str, "_")[[1]]
    if (length(date_parts) > 1) {
      date_code <- date_parts[2] # Extract the date part after the field code
    } else {
      date_code <- date_str # Use whole string if not formatted as expected
    }
    
    # Set band configuration based on date
    # For 7-1-24 and 7-26-24: Band 1-6 = Red, Green, Blue, NIR, Red edge, Deep blue
    # For other dates: Band 1-6 = Deep blue, Blue, Green, Red, Red edge, NIR
    if (grepl("7-1-24", date_code) || grepl("7-26-24", date_code)) {
      band_config <- "early" # Early configuration
      cat(paste("Using early band configuration: Band 1-6 = Red, Green, Blue, NIR, Red edge, Deep blue\n"), 
          file=log_file, append=TRUE)
    } else {
      band_config <- "late" # Later configuration
      cat(paste("Using late band configuration: Band 1-6 = Deep blue, Blue, Green, Red, Red edge, NIR\n"), 
          file=log_file, append=TRUE)
    }
    
    # Process the files
    message("Processing satellite files...")
    results <- process_satellite_plots(satellite_path, shapefile_path, output_folder, band_config)
    
    # Log success
    success_msg <- paste("Successfully processed satellite folder:", date_folder_name)
    cat(paste(success_msg, "\n"), file=log_file, append=TRUE)
    message(success_msg)
    
    # Log completion
    cat(paste("\nProcessing completed at:", Sys.time(), "\n"), file=log_file, append=TRUE)
    cat(sprintf("Total plots: %d\n", nrow(results$raw_data)), file=log_file, append=TRUE)
    cat(sprintf("Total area: %.2f m²\n", sum(results$raw_data$Plot_Area_m2)), file=log_file, append=TRUE)
    message("Processing completed. See log file for details.")
    return(results)
    
  }, error = function(e) {
    # Log any errors that occur during processing
    error_msg <- paste("Error processing satellite folder:", date_folder, "- Error:", conditionMessage(e))
    message(error_msg)
    cat(paste(error_msg, "\n"), file=log_file, append=TRUE)
    return(FALSE)
  })
}

#----------------------------------------------------------------
# SATELLITE PROCESSING FUNCTION
#----------------------------------------------------------------

# Function to process plot data with satellite imagery
process_satellite_plots <- function(satellite_path, shapefile_path, output_folder, band_config) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  # Read raster data
  satellite_raster <- rast(satellite_path)
  
  # Read plot boundaries
  plot_polygons <- st_read(shapefile_path)
  
  # Calculate area for each polygon in square meters
  plot_areas <- st_area(plot_polygons)
  
  # Ensure plot_areas is a numeric vector (not units object)
  plot_areas <- as.numeric(plot_areas)
  
  # Print raster information for debugging
  cat(sprintf("Satellite raster has %d bands\n", nlyr(satellite_raster)))
  cat("Satellite raster band names:", names(satellite_raster), "\n")
  
  # Initialize results dataframe
  id_field <- names(plot_polygons)[grep("id|plot|fid", tolower(names(plot_polygons)))]
  if (length(id_field) == 0) {
    plot_polygons$Plot_ID <- 1:nrow(plot_polygons)
    id_field <- "Plot_ID"
  } else {
    id_field <- id_field[1]
  }
  raw_data <- data.frame(
    Plot_ID = plot_polygons[[id_field]],
    Plot_Area_m2 = plot_areas  # Add area to the raw data
  )
  
  # Log the areas for each plot
  cat("Plot areas (m²):\n")
  for (i in 1:nrow(plot_polygons)) {
    cat(sprintf("Plot %s: %.2f m²\n", plot_polygons[[id_field]][i], plot_areas[i]))
  }
  
  # Create band name mapping based on configuration
  if (band_config == "early") {
    # Early band configuration (7-1-24 and 7-26-24)
    # Band 1-6 = Red, Green, Blue, NIR, Red edge, Deep blue
    satellite_band_mapping <- c(
      "Red" = 1,
      "Green" = 2,
      "Blue" = 3,
      "NIR" = 4,
      "Red_edge" = 5,
      "Deep_blue" = 6
    )
  } else {
    # Late band configuration (other dates)
    # Band 1-6 = Deep blue, Blue, Green, Red, Red edge, NIR
    satellite_band_mapping <- c(
      "Deep_blue" = 1,
      "Blue" = 2,
      "Green" = 3,
      "Red" = 4,
      "Red_edge" = 5,
      "NIR" = 6
    )
  }
  
  # Store pixel-level data for calculating indices
  all_plot_pixels <- list()
  
  # Extract raster statistics for each plot
  cat("Extracting plot statistics...\n")
  
  # Process each plot
  for (i in 1:nrow(plot_polygons)) {
    plot_id <- plot_polygons[[id_field]][i]
    cat(sprintf("Processing plot %d/%d: %s\n", i, nrow(plot_polygons), plot_id))
    
    # Store pixel data for this plot
    all_plot_pixels[[i]] <- list(
      Plot_ID = plot_id,
      Area = plot_areas[i],  # Add the area
      Satellite = list()
    )
    
    # Extract satellite values
    tryCatch({
      satellite_vals <- terra::extract(satellite_raster, plot_polygons[i,])
      # Remove the ID column
      satellite_vals <- satellite_vals[, -1, drop = FALSE]
      
      # Calculate satellite statistics for each band
      for (band_name in names(satellite_band_mapping)) {
        band_idx <- satellite_band_mapping[band_name]
        
        if (band_idx <= ncol(satellite_vals)) {
          band_vals <- satellite_vals[, band_idx]
          
          # Clean values
          band_vals <- band_vals[!is.na(band_vals) & band_vals > 0]
          
          # Store the cleaned band values for index calculation
          all_plot_pixels[[i]]$Satellite[[band_name]] <- band_vals
          
          if (length(band_vals) > 0) {
            # Add statistics to raw_data with meaningful names
            column_prefix <- paste0("Satellite_", band_name)
            raw_data[i, paste0(column_prefix, "_Mean")] <- mean(band_vals)
            raw_data[i, paste0(column_prefix, "_Median")] <- median(band_vals)
            raw_data[i, paste0(column_prefix, "_Sum")] <- sum(band_vals)
            raw_data[i, paste0(column_prefix, "_Sum_per_m2")] <- sum(band_vals) / plot_areas[i]
            # Add pixel count
            raw_data[i, paste0(column_prefix, "_PixelCount")] <- length(band_vals)
          }
        } else {
          cat(sprintf("Warning: Band index %d for %s is out of range for the satellite data\n", 
                      band_idx, band_name))
        }
      }
    }, error = function(e) {
      cat(sprintf("Error processing satellite data for plot %s: %s\n", plot_id, e$message))
    })
  }
  
  # Calculate vegetation indices at the pixel level
  cat("Calculating vegetation indices at pixel level...\n")
  raw_data <- calculate_vegetation_indices_satellite(raw_data, all_plot_pixels, output_folder)
  
  # Process for complete outlier removal with NO NAs
  cat("\nDetecting outliers and creating cleaned dataset with NO NAs...\n")
  cleaned_data_no_nas <- detect_and_completely_remove_outliers(raw_data, output_folder)
  
  return(cleaned_data_no_nas)
}

# Function to calculate vegetation indices at the pixel level for satellite data
calculate_vegetation_indices_satellite <- function(df, all_plot_pixels, output_folder) {
  # Print all column names to debug
  cat("Available column names in dataset:\n")
  cat(paste(names(df), collapse=", "), "\n\n")
  
  # Define platform for satellite
  platform <- "Satellite"
  
  # Function to check if all bands exist for a specific plot's data
  bands_exist_in_plot <- function(plot_data, band_list) {
    all(band_list %in% names(plot_data[[platform]]))
  }
  
  # List of all vegetation indices with their formulas and required bands
  indices <- list(
    # Basic indices
    "NDVI" = list(
      formula = "(NIR - Red) / (NIR + Red)",
      bands = c("NIR", "Red")
    ),
    "GNDVI" = list(
      formula = "(NIR - Green) / (NIR + Green)",
      bands = c("NIR", "Green")
    ),
    "NDRE" = list(
      formula = "(NIR - Red_edge) / (NIR + Red_edge)",
      bands = c("NIR", "Red_edge")
    ),
    "GLI" = list(
      formula = "(2*Green - Red - Blue) / (2*Green + Red + Blue)",
      bands = c("Green", "Red", "Blue")
    ),
    "NGRDI" = list(
      formula = "(Green - Red) / (Green + Red)",
      bands = c("Green", "Red")
    ),
    "SAVI" = list(
      formula = "((NIR - Red) / (NIR + Red + L)) * (1 + L), where L=0.5",
      bands = c("NIR", "Red")
    ),
    "EVI" = list(
      formula = "2.5 * ((NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1))",
      bands = c("NIR", "Red", "Blue")
    ),
    "MSAVI" = list(
      formula = "(2*NIR + 1 - sqrt((2*NIR + 1)^2 - 8*(NIR - Red))) / 2",
      bands = c("NIR", "Red")
    ),
    "NDWI" = list(
      formula = "(Green - NIR) / (Green + NIR)",
      bands = c("Green", "NIR")
    ),
    "SIPI" = list(
      formula = "(NIR - Blue) / (NIR - Red)",
      bands = c("NIR", "Blue", "Red")
    ),
    "MTCI" = list(
      formula = "(NIR - Red_edge) / (Red_edge - Red)",
      bands = c("NIR", "Red_edge", "Red")
    ),
    "CIgreen" = list(
      formula = "(NIR / Green) - 1",
      bands = c("NIR", "Green")
    ),
    "CIrededge" = list(
      formula = "(NIR / Red_edge) - 1",
      bands = c("NIR", "Red_edge")
    ),
    "ARVI" = list(
      formula = "(NIR - (2*Red - Blue)) / (NIR + (2*Red - Blue))",
      bands = c("NIR", "Red", "Blue")
    ),
    "VARI" = list(
      formula = "(Green - Red) / (Green + Red - Blue)",
      bands = c("Green", "Red", "Blue")
    ),
    "OSAVI" = list(
      formula = "(NIR - Red) / (NIR + Red + 0.16)",
      bands = c("NIR", "Red")
    ),
    "TGI" = list(
      formula = "0.5*((670-480)*(Green-Blue)-(550-480)*(Red-Blue))",
      bands = c("Red", "Green", "Blue")
    ),
    "ExG" = list(
      formula = "2*Green - Red - Blue",
      bands = c("Green", "Red", "Blue")
    ),
    "RGRI" = list(
      formula = "Red / Green",
      bands = c("Red", "Green")
    ),
    "TCARI" = list(
      formula = "3 * ((Red_edge - Red) - 0.2 * (Red_edge - Green)) * (Red_edge / Red)",
      bands = c("Red_edge", "Red", "Green")
    )
  )
  
  # Create a summary table to track which indices were calculated
  index_summary <- data.frame(
    Index = character(),
    Status = character(),
    stringsAsFactors = FALSE
  )
  
  # Calculate each index
  for (index_name in names(indices)) {
    index_info <- indices[[index_name]]
    
    # Create a unique name for the index
    index_column_name <- sprintf("%s_%s", index_name, platform)
    cat(sprintf("Calculating %s at pixel level...\n", index_column_name))
    
    # Process each plot
    for (plot_idx in 1:length(all_plot_pixels)) {
      plot_data <- all_plot_pixels[[plot_idx]]
      
      # Check if all required bands exist for this plot
      if (!bands_exist_in_plot(plot_data, index_info$bands)) {
        cat(sprintf("  Skipping plot %s - missing required bands\n", plot_data$Plot_ID))
        next
      }
      
      # Extract band values for this plot
      band_values <- list()
      for (band in index_info$bands) {
        band_values[[band]] <- plot_data[[platform]][[band]]
      }
      
      # Check that all bands have the same number of pixels
      pixel_counts <- sapply(band_values, length)
      if (length(unique(pixel_counts)) != 1) {
        cat(sprintf("  Warning: Plot %s has inconsistent pixel counts across bands\n", plot_data$Plot_ID))
        # Use the minimum number of pixels for all bands
        min_pixels <- min(pixel_counts)
        for (band in index_info$bands) {
          band_values[[band]] <- band_values[[band]][1:min_pixels]
        }
      }
      
      # Calculate index for each pixel
      index_pixels <- NULL
      tryCatch({
        if (index_name == "NDVI") {
          index_pixels <- mapply(function(nir, red) {
            ifelse(!is.na(nir) & !is.na(red) & (nir + red) != 0,
                   (nir - red) / (nir + red), 
                   NA)
          }, band_values$NIR, band_values$Red)
        } else if (index_name == "GNDVI") {
          index_pixels <- mapply(function(nir, green) {
            ifelse(!is.na(nir) & !is.na(green) & (nir + green) != 0,
                   (nir - green) / (nir + green), 
                   NA)
          }, band_values$NIR, band_values$Green)
        } else if (index_name == "NDRE") {
          index_pixels <- mapply(function(nir, rededge) {
            ifelse(!is.na(nir) & !is.na(rededge) & (nir + rededge) != 0,
                   (nir - rededge) / (nir + rededge), 
                   NA)
          }, band_values$NIR, band_values$Red_edge)
        } else if (index_name == "GLI") {
          index_pixels <- mapply(function(green, red, blue) {
            ifelse(!is.na(green) & !is.na(red) & !is.na(blue) & (2*green + red + blue) != 0,
                   (2*green - red - blue) / (2*green + red + blue), 
                   NA)
          }, band_values$Green, band_values$Red, band_values$Blue)
        } else if (index_name == "NGRDI") {
          index_pixels <- mapply(function(green, red) {
            ifelse(!is.na(green) & !is.na(red) & (green + red) != 0,
                   (green - red) / (green + red), 
                   NA)
          }, band_values$Green, band_values$Red)
        } else if (index_name == "SAVI") {
          L <- 0.5
          index_pixels <- mapply(function(nir, red) {
            ifelse(!is.na(nir) & !is.na(red) & (nir + red + L) != 0,
                   ((nir - red) / (nir + red + L)) * (1 + L), 
                   NA)
          }, band_values$NIR, band_values$Red)
        } else if (index_name == "EVI") {
          index_pixels <- mapply(function(nir, red, blue) {
            ifelse(!is.na(nir) & !is.na(red) & !is.na(blue) & (nir + 6*red - 7.5*blue + 1) != 0,
                   2.5 * ((nir - red) / (nir + 6*red - 7.5*blue + 1)),
                   NA)
          }, band_values$NIR, band_values$Red, band_values$Blue)
        } else if (index_name == "MSAVI") {
          index_pixels <- mapply(function(nir, red) {
            ifelse(!is.na(nir) & !is.na(red),
                   (2*nir + 1 - sqrt((2*nir + 1)^2 - 8*(nir - red))) / 2,
                   NA)
          }, band_values$NIR, band_values$Red)
        } else if (index_name == "NDWI") {
          index_pixels <- mapply(function(green, nir) {
            ifelse(!is.na(green) & !is.na(nir) & (green + nir) != 0,
                   (green - nir) / (green + nir), 
                   NA)
          }, band_values$Green, band_values$NIR)
        } else if (index_name == "SIPI") {
          index_pixels <- mapply(function(nir, blue, red) {
            ifelse(!is.na(nir) & !is.na(blue) & !is.na(red) & (nir - red) != 0,
                   (nir - blue) / (nir - red), 
                   NA)
          }, band_values$NIR, band_values$Blue, band_values$Red)
        } else if (index_name == "MTCI") {
          index_pixels <- mapply(function(nir, rededge, red) {
            ifelse(!is.na(nir) & !is.na(rededge) & !is.na(red) & (rededge - red) != 0,
                   (nir - rededge) / (rededge - red), 
                   NA)
          }, band_values$NIR, band_values$Red_edge, band_values$Red)
        } else if (index_name == "CIgreen") {
          index_pixels <- mapply(function(nir, green) {
            ifelse(!is.na(nir) & !is.na(green) & green != 0,
                   (nir / green) - 1, 
                   NA)
          }, band_values$NIR, band_values$Green)
        } else if (index_name == "CIrededge") {
          index_pixels <- mapply(function(nir, rededge) {
            ifelse(!is.na(nir) & !is.na(rededge) & rededge != 0,
                   (nir / rededge) - 1, 
                   NA)
          }, band_values$NIR, band_values$Red_edge)
        } else if (index_name == "ARVI") {
          index_pixels <- mapply(function(nir, red, blue) {
            rb <- 2*red - blue
            ifelse(!is.na(nir) & !is.na(red) & !is.na(blue) & (nir + rb) != 0,
                   (nir - rb) / (nir + rb), 
                   NA)
          }, band_values$NIR, band_values$Red, band_values$Blue)
        } else if (index_name == "VARI") {
          index_pixels <- mapply(function(green, red, blue) {
            ifelse(!is.na(green) & !is.na(red) & !is.na(blue) & (green + red - blue) != 0,
                   (green - red) / (green + red - blue), 
                   NA)
          }, band_values$Green, band_values$Red, band_values$Blue)
        } else if (index_name == "OSAVI") {
          index_pixels <- mapply(function(nir, red) {
            ifelse(!is.na(nir) & !is.na(red) & (nir + red + 0.16) != 0,
                   (nir - red) / (nir + red + 0.16), 
                   NA)
          }, band_values$NIR, band_values$Red)
        } else if (index_name == "TGI") {
          index_pixels <- mapply(function(red, green, blue) {
            ifelse(!is.na(red) & !is.na(green) & !is.na(blue),
                   0.5*((670-480)*(green-blue)-(550-480)*(red-blue)), 
                   NA)
          }, band_values$Red, band_values$Green, band_values$Blue)
        } else if (index_name == "ExG") {
          index_pixels <- mapply(function(green, red, blue) {
            ifelse(!is.na(green) & !is.na(red) & !is.na(blue),
                   2*green - red - blue, 
                   NA)
          }, band_values$Green, band_values$Red, band_values$Blue)
        } else if (index_name == "RGRI") {
          index_pixels <- mapply(function(red, green) {
            ifelse(!is.na(red) & !is.na(green) & green != 0,
                   red / green, 
                   NA)
          }, band_values$Red, band_values$Green)
        } else if (index_name == "TCARI") {
          index_pixels <- mapply(function(rededge, red, green) {
            ifelse(!is.na(rededge) & !is.na(red) & !is.na(green) & red != 0,
                   3 * ((rededge - red) - 0.2 * (rededge - green) * (rededge / red)), 
                   NA)
          }, band_values$Red_edge, band_values$Red, band_values$Green)
        }
        
        # Clean up index_pixels
        index_pixels <- index_pixels[!is.na(index_pixels) & is.finite(index_pixels)]
        
        # Calculate statistics from pixel-level indices
        if (length(index_pixels) > 0) {
          row_idx <- plot_idx  # The row index is the same as the plot index
          
          # Calculate true statistics from pixel-level indices
          df[row_idx, paste0(index_column_name, "_Mean")] <- mean(index_pixels)
          df[row_idx, paste0(index_column_name, "_Median")] <- median(index_pixels)
          df[row_idx, paste0(index_column_name, "_Sum")] <- sum(index_pixels)
          df[row_idx, paste0(index_column_name, "_Sum_per_m2")] <- sum(index_pixels) / all_plot_pixels[[plot_idx]]$Area
          # Add pixel count
          df[row_idx, paste0(index_column_name, "_PixelCount")] <- length(index_pixels)
          
          # Add success to index summary
          if (plot_idx == 1) {  # Only add to summary once per index
            index_summary <- rbind(index_summary, data.frame(
              Index = index_name,
              Status = "Successfully calculated at pixel level",
              stringsAsFactors = FALSE
            ))
          }
        }
      }, error = function(e) {
        cat(sprintf("  Error calculating %s for plot %s: %s\n", 
                    index_name, plot_data$Plot_ID, e$message))
        
        # Add error to index summary only once per index
        if (plot_idx == 1) {
          index_summary <- rbind(index_summary, data.frame(
            Index = index_name,
            Status = sprintf("Error: %s", e$message),
            stringsAsFactors = FALSE
          ))
        }
      })
    }
  }
  
  # Print summary of calculated indices
  cat("\n===== VEGETATION INDICES CALCULATION SUMMARY =====\n")
  summary_table <- table(index_summary$Status)
  for (status in names(summary_table)) {
    cat(sprintf("%s: %d\n", status, summary_table[status]))
  }
  
  cat("\nSuccessfully calculated indices:\n")
  success_indices <- index_summary[grepl("Successfully", index_summary$Status), "Index"]
  cat(paste("  ", success_indices, collapse = "\n  "), "\n")
  
  # Save the index summary to CSV
  index_summary_path <- file.path(output_folder, "idx.csv")
  tryCatch({
    write.csv(index_summary, index_summary_path, row.names = FALSE)
    cat(sprintf("Saved vegetation indices summary to %s\n", index_summary_path))
  }, error = function(e) {
    cat(sprintf("Warning: Could not save vegetation indices summary: %s\n", e$message))
  })
  
  return(df)
}

# Function to detect outliers and completely remove them (no NAs)
detect_and_completely_remove_outliers <- function(df, output_folder) {
  cat(sprintf("Analyzing %d plot records with %d variables\n", nrow(df), ncol(df)))
  
  # Save the raw data with original statistics first
  raw_path <- file.path(output_folder, "raw.csv")
  tryCatch({
    write.csv(df, raw_path, row.names = FALSE)
    cat(sprintf("Saved raw dataset to %s\n", raw_path))
  }, error = function(e) {
    cat(sprintf("Warning: Could not save raw dataset: %s\n", e$message))
  })
  
  # Identify indices and band statistics
  # Look for vegetation indices and their stats (Mean, Median, Sum)
  index_pattern <- "^[A-Za-z0-9]+_Satellite(_Mean|_Median|_Sum)?$"
  index_cols <- grep(index_pattern, names(df), value = TRUE)
  
  # Get band statistics columns - looking for Mean, Median, Sum
  band_stats_pattern <- "_(Mean|Median|Sum|PixelCount)$"  # Added PixelCount to the pattern
  band_stats_cols <- names(df)[grep(band_stats_pattern, names(df))]
  
  # Filter out index statistic columns that are already in index_cols
  band_stats_cols <- setdiff(band_stats_cols, index_cols)
  
  # All columns to analyze
  analyze_cols <- c(index_cols, band_stats_cols)
  cat(sprintf("Found %d vegetation indices and their stats, and %d band statistics to analyze for outliers\n", 
              length(index_cols), length(band_stats_cols)))
  
  # Initialize outlier tracker
  outlier_tracker <- data.frame(Plot_ID = df$Plot_ID)
  outlier_count <- rep(0, nrow(df))
  
  # Track number of outliers per column
  col_outliers <- data.frame(
    column = character(),
    n_outliers = integer(),
    outlier_percent = numeric(),
    stringsAsFactors = FALSE
  )
  
  # First pass - identify outliers in each column
  for (col in analyze_cols) {
    if (!col %in% names(df)) {
      cat(sprintf("Warning: Column %s not found in dataset. Skipping.\n", col))
      next
    }
    
    # Skip columns with all NA values
    if (all(is.na(df[[col]]))) {
      cat(sprintf("Warning: Column %s contains all NA values. Skipping.\n", col))
      next
    }
    
    # Get column data, dropping NA values
    col_data <- df[[col]][!is.na(df[[col]])]
    
    if (length(col_data) == 0) {
      cat(sprintf("Warning: No valid data in column %s after dropping NAs. Skipping.\n", col))
      next
    }
    
    # Calculate z-scores (handle single-value columns)
    if (length(unique(col_data)) == 1) {
      z_scores <- rep(0, length(col_data))  # No variation, so no outliers
    } else {
      z_scores <- abs(scale(col_data))
    }
    
    # Identify outliers (z-score > 3)
    outlier_indices <- which(z_scores > 3)
    
    # Track which plots have outliers in this column
    outlier_tracker[[paste0(col, "_outlier")]] <- FALSE
    
    if (length(outlier_indices) > 0) {
      original_indices <- which(!is.na(df[[col]]))
      plot_indices_with_outliers <- original_indices[outlier_indices]
      
      # Update outlier count for these plots
      outlier_count[plot_indices_with_outliers] <- outlier_count[plot_indices_with_outliers] + 1
      
      if (length(plot_indices_with_outliers) > 0) {
        outlier_tracker[plot_indices_with_outliers, paste0(col, "_outlier")] <- TRUE
      }
      
      # Add to column outlier tracking
      col_outliers <- rbind(col_outliers, data.frame(
        column = col,
        n_outliers = length(outlier_indices),
        outlier_percent = 100 * length(outlier_indices) / length(col_data),
        stringsAsFactors = FALSE
      ))
    } else {
      # No outliers in this column
      col_outliers <- rbind(col_outliers, data.frame(
        column = col,
        n_outliers = 0,
        outlier_percent = 0,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Add outlier count to tracker
  outlier_tracker$outlier_count <- outlier_count
  outlier_tracker$Plot_Area_m2 <- df$Plot_Area_m2  # Add plot area to the tracker
  
  # Save outlier tracker
  outlier_tracker_path <- file.path(output_folder, "out.csv")
  tryCatch({
    write.csv(outlier_tracker, outlier_tracker_path, row.names = FALSE)
    cat(sprintf("Saved outlier tracker to %s\n", outlier_tracker_path))
  }, error = function(e) {
    cat(sprintf("Warning: Could not save outlier tracker: %s\n", e$message))
  })
  
  # Create data frame of just the problematic plots (those with any outliers)
  problem_plots <- outlier_tracker$Plot_ID[outlier_tracker$outlier_count > 0]
  cat(sprintf("\nFound %d plots with outliers (%.1f%% of total)\n", 
              length(problem_plots), 100 * length(problem_plots) / nrow(df)))
  
  # Create cleaned dataset by removing problem plots
  df_cleaned <- df[!df$Plot_ID %in% problem_plots, ]
  cat(sprintf("Successfully created cleaned dataset with %d plots (removed %d outlier plots)\n", 
              nrow(df_cleaned), length(problem_plots)))
  
  # Save cleaned dataset (no NAs, completely removed outlier plots)
  cleaned_path <- file.path(output_folder, "cln.csv")
  tryCatch({
    write.csv(df_cleaned, cleaned_path, row.names = FALSE)
    cat(sprintf("Saved cleaned dataset to %s\n", cleaned_path))
  }, error = function(e) {
    cat(sprintf("Warning: Could not save cleaned dataset: %s\n", e$message))
  })
  
  # Calculate and display overall statistics
  total_values <- sum(sapply(df[analyze_cols], function(x) sum(!is.na(x))))
  total_outliers <- sum(col_outliers$n_outliers)
  overall_outlier_percent <- 100 * total_outliers / total_values
  
  # Calculate area statistics
  total_area <- sum(df$Plot_Area_m2)
  cleaned_area <- sum(df_cleaned$Plot_Area_m2)
  area_percent_remaining <- 100 * cleaned_area / total_area
  
  # Create summary report for saving
  summary_report <- data.frame(
    metric = c(
      "Total_data_points_analyzed", 
      "Total_outliers_detected", 
      "Overall_outlier_percentage",
      "Plots_with_outliers",
      "Plots_with_outliers_percent",
      "Plots_in_cleaned_dataset",
      "Total_area_m2",
      "Cleaned_area_m2",
      "Area_percent_remaining"
    ),
    value = c(
      total_values,
      total_outliers,
      overall_outlier_percent,
      length(problem_plots),
      100 * length(problem_plots) / nrow(df),
      nrow(df_cleaned),
      total_area,
      cleaned_area,
      area_percent_remaining
    )
  )
  
  # Save the outlier summary report
  summary_path <- file.path(output_folder, "sum.csv")
  tryCatch({
    write.csv(summary_report, summary_path, row.names = FALSE)
    cat(sprintf("Saved outlier summary report to %s\n", summary_path))
    
    # Print top variables with highest outlier percentages (with safety checks)
    cat("\nTop variables with highest outlier percentages:\n")
    if (nrow(col_outliers) > 0) {
      # Sort outlier summary by outlier_percent in descending order
      sorted_summary <- col_outliers[order(col_outliers$outlier_percent, decreasing = TRUE), ]
      
      # Get the top 10 or fewer if there are less than 10 rows
      top_n <- min(10, nrow(sorted_summary))
      if (top_n > 0) {
        top_outliers <- sorted_summary[1:top_n, ]
        
        for (i in 1:nrow(top_outliers)) {
          cat(sprintf("  %s: %.2f%% (%d outliers)\n", 
                      top_outliers$column[i], top_outliers$outlier_percent[i], top_outliers$n_outliers[i]))
        }
        
        # Add top outliers to the summary report
        top_outliers_report <- data.frame(
          metric = paste0("Top_", 1:top_n),
          value = top_outliers$column,
          percent = top_outliers$outlier_percent,
          count = top_outliers$n_outliers
        )
        
        # Append to the existing CSV
        tryCatch({
          write.table(top_outliers_report, summary_path, append = TRUE, 
                      sep = ",", col.names = TRUE, row.names = FALSE)
        }, error = function(e) {
          cat(sprintf("Warning: Could not append top outliers report: %s\n", e$message))
        })
      } else {
        cat("  No variables with outliers found.\n")
      }
    } else {
      cat("  No outliers detected in any variables.\n")
    }
  }, error = function(e) {
    cat(sprintf("Warning: Could not save outlier summary report: %s\n", e$message))
  })
  
  # Print statistics to console
  cat("\n===== OUTLIER DETECTION SUMMARY =====\n")
  cat(sprintf("Total data points analyzed: %d\n", total_values))
  cat(sprintf("Total outliers detected: %d\n", total_outliers))
  cat(sprintf("Overall outlier percentage: %.2f%%\n", overall_outlier_percent))
  cat(sprintf("Plots with at least one outlier: %d out of %d (%.2f%%)\n", 
              length(problem_plots), nrow(df), 100*length(problem_plots)/nrow(df)))
  cat(sprintf("Plots in cleaned dataset: %d\n", nrow(df_cleaned)))
  cat(sprintf("Total area: %.2f m²\n", total_area))
  cat(sprintf("Cleaned area: %.2f m² (%.2f%% of total)\n", cleaned_area, area_percent_remaining))
  
  return(list(
    raw_data = df,
    cleaned_data = df_cleaned,
    outlier_summary = col_outliers,
    summary_report = summary_report,
    outlier_tracker = outlier_tracker
  ))
}

#----------------------------------------------------------------
# FUNCTIONS TO RETRIEVE SATELLITE FOLDERS
#----------------------------------------------------------------

# Function to get all Satellite folders
get_all_satellite_folders <- function(base_dir) {
  folder_list <- list()
  
  for (field in c("EF", "SF")) {
    field_path <- file.path(base_dir, field)
    if (!dir.exists(field_path)) {
      message(paste("Field directory not found:", field_path))
      next
    }
    
    satellite_path <- file.path(field_path, "Satellite")
    if (!dir.exists(satellite_path)) {
      message(paste("Satellite directory not found for field:", field))
      next
    }
    
    date_folders <- list.dirs(satellite_path, full.names = FALSE, recursive = FALSE)
    
    if (length(date_folders) > 0) {
      for (folder in date_folders) {
        folder_list <- c(folder_list, list(list(field = field, folder = folder)))
      }
    }
  }
  
  return(folder_list)
}

#----------------------------------------------------------------
# MAIN BATCH PROCESSING FUNCTION FOR SATELLITE DATA
#----------------------------------------------------------------

# Function to process all satellite folders sequentially
process_all_satellite_folders_sequentially <- function(base_dir = getwd()) {
  # Create main log file
  master_log <- file.path(base_dir, "satellite_batch_processing_master_log.txt")
  cat(paste("Satellite Batch Processing started at:", Sys.time(), "\n"), file = master_log)
  
  # Get all Satellite folders to process
  folders <- get_all_satellite_folders(base_dir)
  
  if (length(folders) == 0) {
    message("No Satellite folders found to process.")
    cat("No Satellite folders found to process.\n", file = master_log, append = TRUE)
    return()
  }
  
  cat(sprintf("Found %d Satellite folders to process.\n", length(folders)), file = master_log, append = TRUE)
  
  # Create a summary table
  summary_table <- data.frame(
    Field = character(),
    Folder = character(),
    Status = character(),
    ProcessingTime = numeric(),
    PlotCount = integer(),
    TotalArea_m2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each folder sequentially
  for (i in seq_along(folders)) {
    field <- folders[[i]]$field
    folder <- folders[[i]]$folder
    
    message(sprintf("\nProcessing satellite folder %d of %d: %s/%s", i, length(folders), field, folder))
    cat(sprintf("\nProcessing satellite folder %d of %d: %s/%s at %s\n", 
                i, length(folders), field, folder, Sys.time()), 
        file = master_log, append = TRUE)
    
    # Track processing time
    start_time <- Sys.time()
    
    # Process the folder
    result <- tryCatch({
      output <- process_single_satellite_folder(base_dir, field, folder)
      # Check if result is a list containing data (successful run with returned data)
      if (is.list(output) && "raw_data" %in% names(output)) {
        list(success = TRUE, plot_count = nrow(output$raw_data), total_area = sum(output$raw_data$Plot_Area_m2))
      } else {
        list(success = output, plot_count = NA, total_area = NA)
      }
    }, error = function(e) {
      error_msg <- sprintf("Error in main processing loop: %s", conditionMessage(e))
      message(error_msg)
      cat(paste(error_msg, "\n"), file = master_log, append = TRUE)
      return(list(success = FALSE, plot_count = NA, total_area = NA))
    })
    
    # Calculate processing time
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
    
    # Update summary table
    summary_table <- rbind(summary_table, data.frame(
      Field = field,
      Folder = folder,
      Status = ifelse(result$success, "Success", "Failed"),
      ProcessingTime = processing_time,
      PlotCount = ifelse(is.na(result$plot_count), 0, result$plot_count),
      TotalArea_m2 = ifelse(is.na(result$total_area), 0, result$total_area),
      stringsAsFactors = FALSE
    ))
    
    # Log result
    result_msg <- sprintf("Satellite folder %s/%s processed %s in %.2f minutes", 
                          field, folder, 
                          ifelse(result$success, "successfully", "with errors"),
                          processing_time)
    message(result_msg)
    cat(paste(result_msg, "\n"), file = master_log, append = TRUE)
    
    # Save current progress to summary CSV
    summary_path <- file.path(base_dir, "satellite_batch_processing_summary.csv")
    write.csv(summary_table, summary_path, row.names = FALSE)
    
    # Run garbage collection to free memory before next folder
    gc()
  }
  
  # Print final summary
  success_count <- sum(summary_table$Status == "Success")
  failed_count <- sum(summary_table$Status == "Failed")
  total_time <- sum(summary_table$ProcessingTime)
  total_plots <- sum(summary_table$PlotCount, na.rm = TRUE)
  total_area <- sum(summary_table$TotalArea_m2, na.rm = TRUE)
  
  summary_msg <- sprintf("\nSatellite batch processing completed at %s\n", Sys.time())
  summary_msg <- paste0(summary_msg, sprintf("Successfully processed: %d folders\n", success_count))
  summary_msg <- paste0(summary_msg, sprintf("Failed: %d folders\n", failed_count))
  summary_msg <- paste0(summary_msg, sprintf("Total processing time: %.2f minutes\n", total_time))
  summary_msg <- paste0(summary_msg, sprintf("Total plots processed: %d\n", total_plots))
  summary_msg <- paste0(summary_msg, sprintf("Total area processed: %.2f m²\n", total_area))
  
  message(summary_msg)
  cat(summary_msg, file = master_log, append = TRUE)
  
  # List failed folders if any
  if (failed_count > 0) {
    failed_folders <- summary_table[summary_table$Status == "Failed", c("Field", "Folder")]
    failed_list <- paste(sprintf("%s/%s", failed_folders$Field, failed_folders$Folder), collapse = "\n  ")
    failed_msg <- sprintf("\nFailed satellite folders:\n  %s\n", failed_list)
    message(failed_msg)
    cat(failed_msg, file = master_log, append = TRUE)
  }
}

#----------------------------------------------------------------
# RUN THE BATCH PROCESSING FOR SATELLITE DATA
#----------------------------------------------------------------

# Get the base directory
base_dir <- getwd()
cat("Base directory:", base_dir, "\n")

# Run the batch processing for satellite data
cat("\n========== STARTING SATELLITE BATCH PROCESSING ==========\n")
process_all_satellite_folders_sequentially(base_dir)
cat("\n========== SATELLITE BATCH PROCESSING COMPLETE ==========\n")