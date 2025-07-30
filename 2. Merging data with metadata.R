# Adaptive Data Processing Script
# This script merges data files with metadata in a flexible, configurable way

# Load required libraries
required_packages <- c("readr", "dplyr", "tidyr", "readxl", "stringr", "yaml", "here")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Configuration class to handle all settings
DataProcessorConfig <- list(
  # Directory settings - automatically detects project structure
  # When script is run from Codes folder, base_dir points to parent
  base_dir = here::here(".."),  # Go up one level from Codes folder to project root
  split_data_dir = NULL,  # Will default to ../Analysis/2. Split_Data
  metadata_dir = NULL,    # Will default to ../Data/Orthomosaics+shapefile/Metadata
  output_dir = NULL,      # Will default to ../Analysis/4. Merged_data
  
  # File pattern settings
  data_file_pattern = "_data\\.csv$",
  metadata_file_name = "Soybean_UAV_Satellite_metadata.xlsx",
  metadata_sheet_name = "Ground_truth",
  
  # Column settings
  required_data_columns = c("Plot_ID"),
  required_metadata_columns = c("Plot_ID"),
  metadata_columns_to_keep = NULL, # NULL means keep all
  
  # Joining settings
  plot_id_column = "Plot_ID",
  location_column = "Location",
  
  # File naming conventions
  filename_separator = "_",
  location_position = 1,
  platform_position = 2,
  
  # Processing options
  case_sensitive_ids = FALSE,
  auto_detect_prefixes = TRUE,
  verbose = TRUE,
  save_log = TRUE
)

#' Load configuration from file or use defaults
#' @param config_file Path to YAML configuration file (optional)
#' @return Updated configuration list
load_config <- function(config_file = NULL) {
  config <- DataProcessorConfig
  
  if (!is.null(config_file) && file.exists(config_file)) {
    user_config <- yaml::read_yaml(config_file)
    config <- modifyList(config, user_config)
    if (config$verbose) cat("Configuration loaded from:", config_file, "\n")
  }
  
  # Set default directories if not specified
  # Assumes script is run from Codes folder, so base_dir is parent directory
  if (is.null(config$split_data_dir)) {
    config$split_data_dir <- file.path(config$base_dir, "Analysis", "2. Split_Data")
  }
  if (is.null(config$metadata_dir)) {
    config$metadata_dir <- file.path(config$base_dir, "Data", "Orthomosaics+shapefile", "Metadata")
  }
  if (is.null(config$output_dir)) {
    config$output_dir <- file.path(config$base_dir, "Analysis", "4. Merged_data")
  }
  
  return(config)
}

#' Validate configuration and directories
#' @param config Configuration list
#' @return TRUE if valid, stops execution if invalid
validate_config <- function(config) {
  # Check if directories exist
  dirs_to_check <- c("split_data_dir", "metadata_dir")
  for (dir_name in dirs_to_check) {
    dir_path <- config[[dir_name]]
    if (!dir.exists(dir_path)) {
      stop(paste("Directory does not exist:", dir_path, 
                 "\nPlease check your configuration or create the directory."))
    }
  }
  
  # Create output directory if it doesn't exist
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Check if metadata file exists
  metadata_path <- file.path(config$metadata_dir, config$metadata_file_name)
  if (!file.exists(metadata_path)) {
    stop(paste("Metadata file not found:", metadata_path))
  }
  
  return(TRUE)
}

#' Extract location and platform information from filename
#' @param filename The filename to parse
#' @param config Configuration list
#' @return List with location and platform, or NULL if parsing fails
extract_file_info <- function(filename, config) {
  # Remove file extension
  base_name <- tools::file_path_sans_ext(basename(filename))
  
  # Split by separator
  parts <- str_split(base_name, config$filename_separator, simplify = TRUE)
  
  if (length(parts) >= max(config$location_position, config$platform_position)) {
    location <- parts[config$location_position]
    platform <- parts[config$platform_position]
    
    # Handle case sensitivity
    if (!config$case_sensitive_ids) {
      location <- toupper(location)
      platform <- toupper(platform)
    }
    
    return(list(location = location, platform = platform))
  }
  
  return(NULL)
}

#' Load and validate metadata
#' @param config Configuration list
#' @return Processed metadata dataframe
load_metadata <- function(config) {
  metadata_path <- file.path(config$metadata_dir, config$metadata_file_name)
  
  if (config$verbose) cat("Loading metadata from:", metadata_path, "\n")
  
  # Try to read metadata with error handling
  tryCatch({
    if (grepl("\\.xlsx?$", config$metadata_file_name, ignore.case = TRUE)) {
      metadata <- read_excel(metadata_path, sheet = config$metadata_sheet_name)
    } else if (grepl("\\.csv$", config$metadata_file_name, ignore.case = TRUE)) {
      metadata <- read_csv(metadata_path, show_col_types = FALSE)
    } else {
      stop("Unsupported metadata file format. Use .xlsx, .xls, or .csv")
    }
  }, error = function(e) {
    stop(paste("Error reading metadata file:", e$message))
  })
  
  # Validate required columns
  missing_cols <- setdiff(config$required_metadata_columns, names(metadata))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in metadata:", paste(missing_cols, collapse = ", ")))
  }
  
  # Select columns to keep
  if (!is.null(config$metadata_columns_to_keep)) {
    available_cols <- intersect(config$metadata_columns_to_keep, names(metadata))
    metadata <- metadata %>% select(all_of(available_cols))
  }
  
  # Ensure Plot_ID is character for consistent joining
  metadata[[config$plot_id_column]] <- as.character(metadata[[config$plot_id_column]])
  
  # Handle case sensitivity
  if (!config$case_sensitive_ids) {
    if (config$location_column %in% names(metadata)) {
      metadata[[config$location_column]] <- toupper(metadata[[config$location_column]])
    }
  }
  
  if (config$verbose) {
    cat("Metadata loaded successfully. Rows:", nrow(metadata), "Columns:", ncol(metadata), "\n")
  }
  
  return(metadata)
}

#' Find potential ID matching strategies
#' @param data_ids Vector of Plot_IDs from data file
#' @param metadata_ids Vector of Plot_IDs from metadata
#' @param location Current location being processed
#' @param config Configuration list
#' @return List with matching strategy and transformed IDs
find_id_matching_strategy <- function(data_ids, metadata_ids, location, config) {
  # Direct match
  common_ids <- intersect(data_ids, metadata_ids)
  if (length(common_ids) > 0) {
    return(list(strategy = "direct", join_ids = data_ids, matches = length(common_ids)))
  }
  
  if (!config$auto_detect_prefixes) {
    return(list(strategy = "none", join_ids = data_ids, matches = 0))
  }
  
  # Try various prefix/suffix patterns
  patterns_to_try <- c(
    paste0(location, "_"),
    paste0(location, "."),
    paste0(tolower(location), "_"),
    paste0(tolower(location), "."),
    paste0(location, "-"),
    "_", "."
  )
  
  for (pattern in patterns_to_try) {
    # Try adding prefix to data IDs
    prefixed_data_ids <- paste0(pattern, data_ids)
    common_ids <- intersect(prefixed_data_ids, metadata_ids)
    if (length(common_ids) > 0) {
      return(list(strategy = paste("add_prefix:", pattern), 
                  join_ids = prefixed_data_ids, 
                  matches = length(common_ids)))
    }
    
    # Try removing prefix from data IDs
    if (any(grepl(paste0("^", pattern), data_ids))) {
      stripped_data_ids <- gsub(paste0("^", pattern), "", data_ids)
      common_ids <- intersect(stripped_data_ids, metadata_ids)
      if (length(common_ids) > 0) {
        return(list(strategy = paste("remove_prefix:", pattern), 
                    join_ids = stripped_data_ids, 
                    matches = length(common_ids)))
      }
    }
    
    # Try removing prefix from metadata IDs (conceptually)
    if (any(grepl(paste0("^", pattern), metadata_ids))) {
      stripped_metadata_ids <- gsub(paste0("^", pattern), "", metadata_ids)
      common_ids <- intersect(data_ids, stripped_metadata_ids)
      if (length(common_ids) > 0) {
        return(list(strategy = paste("metadata_has_prefix:", pattern), 
                    join_ids = data_ids, 
                    matches = length(common_ids),
                    metadata_transform = TRUE))
      }
    }
  }
  
  return(list(strategy = "none", join_ids = data_ids, matches = 0))
}

#' Process a single data file
#' @param file_path Path to the data file
#' @param metadata Complete metadata dataframe
#' @param config Configuration list
#' @return Processed and merged dataframe, or NULL if processing fails
process_data_file <- function(file_path, metadata, config) {
  filename <- basename(file_path)
  
  # Extract file information
  file_info <- extract_file_info(file_path, config)
  if (is.null(file_info)) {
    if (config$verbose) {
      cat("WARNING: Could not parse location and platform from filename:", filename, "\n")
    }
    return(NULL)
  }
  
  location <- file_info$location
  platform <- file_info$platform
  
  if (config$verbose) cat("\nProcessing", location, platform, "data...\n")
  
  # Read data file with error handling
  tryCatch({
    data <- read_csv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    if (config$verbose) {
      cat("ERROR: Could not read data file:", filename, "-", e$message, "\n")
    }
    return(NULL)
  })
  
  if (is.null(data)) return(NULL)
  
  # Validate required columns
  missing_cols <- setdiff(config$required_data_columns, names(data))
  if (length(missing_cols) > 0) {
    if (config$verbose) {
      cat("ERROR: Missing required columns in data file:", paste(missing_cols, collapse = ", "), "\n")
    }
    return(NULL)
  }
  
  # Ensure Plot_ID is character
  data[[config$plot_id_column]] <- as.character(data[[config$plot_id_column]])
  
  # Filter metadata for current location
  if (config$location_column %in% names(metadata)) {
    location_metadata <- metadata %>% filter(.data[[config$location_column]] == location)
  } else {
    location_metadata <- metadata
  }
  
  if (config$verbose) {
    cat("  Metadata rows for location", location, ":", nrow(location_metadata), "\n")
  }
  
  if (nrow(location_metadata) == 0) {
    if (config$verbose) {
      cat("  WARNING: No metadata found for location:", location, "\n")
    }
    return(NULL)
  }
  
  # Find matching strategy
  matching_result <- find_id_matching_strategy(
    data[[config$plot_id_column]], 
    location_metadata[[config$plot_id_column]], 
    location, 
    config
  )
  
  if (config$verbose) {
    cat("  Matching strategy:", matching_result$strategy, "\n")
    cat("  Number of matches found:", matching_result$matches, "\n")
  }
  
  if (matching_result$matches == 0) {
    if (config$verbose) {
      cat("  ERROR: No matching Plot_IDs found! Skipping file.\n")
      cat("  Sample data Plot_IDs:", paste(head(unique(data[[config$plot_id_column]]), 3), collapse = ", "), "\n")
      cat("  Sample metadata Plot_IDs:", paste(head(unique(location_metadata[[config$plot_id_column]]), 3), collapse = ", "), "\n")
    }
    return(NULL)
  }
  
  # Prepare data for joining
  data$Plot_ID_join <- matching_result$join_ids
  
  # Handle metadata transformation if needed
  if (!is.null(matching_result$metadata_transform) && matching_result$metadata_transform) {
    # This is a placeholder for more complex metadata transformations if needed
    join_column <- config$plot_id_column
  } else {
    join_column <- config$plot_id_column
  }
  
  # Perform the merge
  combined_data <- data %>%
    inner_join(location_metadata, by = c("Plot_ID_join" = join_column)) %>%
    select(-Plot_ID_join) %>%
    mutate(Platform = platform, 
           Location = location,
           ProcessedFile = filename)
  
  if (config$verbose) {
    cat("  Combined rows:", nrow(combined_data), "\n")
    cat("  Unique Plot_IDs in combined data:", length(unique(combined_data[[config$plot_id_column]])), "\n")
  }
  
  return(combined_data)
}

#' Main processing function
#' @param config_file Optional path to configuration file
#' @return List with processing results
process_all_data <- function(config_file = NULL) {
  # Load configuration
  config <- load_config(config_file)
  
  # Validate configuration
  validate_config(config)
  
  # Initialize log
  log_entries <- character()
  log_entry <- function(message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    entry <- paste0("[", timestamp, "] ", message)
    log_entries <<- c(log_entries, entry)
    if (config$verbose) cat(entry, "\n")
  }
  
  log_entry("Starting data processing...")
  log_entry(paste("Data directory:", config$split_data_dir))
  log_entry(paste("Metadata directory:", config$metadata_dir))
  log_entry(paste("Output directory:", config$output_dir))
  
  # Load metadata
  metadata <- load_metadata(config)
  log_entry(paste("Loaded metadata with", nrow(metadata), "rows"))
  
  # Get list of data files
  data_files <- list.files(config$split_data_dir, 
                           pattern = config$data_file_pattern, 
                           full.names = TRUE)
  
  log_entry(paste("Found", length(data_files), "data files to process"))
  
  if (length(data_files) == 0) {
    log_entry("No data files found matching the pattern!")
    return(list(success = FALSE, message = "No data files found", log = log_entries))
  }
  
  # Process each file
  successful_files <- character()
  failed_files <- character()
  
  for (file_path in data_files) {
    filename <- basename(file_path)
    log_entry(paste("Processing file:", filename))
    
    combined_data <- process_data_file(file_path, metadata, config)
    
    if (is.null(combined_data)) {
      failed_files <- c(failed_files, filename)
      log_entry(paste("FAILED to process:", filename))
      next
    }
    
    # Generate output filename
    file_info <- extract_file_info(file_path, config)
    output_filename <- paste0(file_info$location, "_", file_info$platform, "_combined.csv")
    output_path <- file.path(config$output_dir, output_filename)
    
    # Save combined data
    tryCatch({
      write_csv(combined_data, output_path)
      successful_files <- c(successful_files, filename)
      log_entry(paste("Successfully saved:", output_filename))
    }, error = function(e) {
      failed_files <- c(failed_files, filename)
      log_entry(paste("ERROR saving file:", e$message))
    })
  }
  
  # Summary
  log_entry("Processing complete!")
  log_entry(paste("Successful files:", length(successful_files)))
  log_entry(paste("Failed files:", length(failed_files)))
  
  # Save log if requested
  if (config$save_log) {
    log_file <- file.path(config$output_dir, paste0("processing_log_", 
                                                    format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
    writeLines(log_entries, log_file)
    cat("Log saved to:", log_file, "\n")
  }
  
  return(list(
    success = length(failed_files) == 0,
    successful_files = successful_files,
    failed_files = failed_files,
    log = log_entries,
    config = config
  ))
}

# Example usage and configuration creation
create_example_config <- function(config_path = "data_processor_config.yaml") {
  example_config <- list(
    # Directory settings for your project structure
    # Script assumes it's run from Codes folder
    base_dir = here::here(".."),  # Parent directory containing Codes, Analysis, Data
    split_data_dir = NULL,  # Will use Analysis/2. Split_Data (auto-detected)
    metadata_dir = NULL,    # Will use Data/Orthomosaics+shapefile/Metadata (auto-detected)
    output_dir = NULL,      # Will use Analysis/4. Merged_data (auto-detected)
    
    # File settings
    data_file_pattern = "_data\\.csv$",
    metadata_file_name = "Soybean_UAV_Satellite_metadata.xlsx",
    metadata_sheet_name = "Ground_truth",
    
    # Columns to keep from metadata (NULL means keep all)
    metadata_columns_to_keep = c("Location", "Loc_ID", "Plot_ID", "Line", "Spacing", 
                                 "Leaf_type", "area", "lai", "lr", "lma", "mxwidth", 
                                 "protein", "oil", "yield"),
    
    # Processing options
    case_sensitive_ids = FALSE,
    auto_detect_prefixes = TRUE,
    verbose = TRUE,
    save_log = TRUE,
    
    # File naming convention
    filename_separator = "_",
    location_position = 1,
    platform_position = 2
  )
  
  yaml::write_yaml(example_config, config_path)
  cat("Example configuration saved to:", config_path, "\n")
  return(config_path)
}

# Run the processing (uncomment to execute)
# IMPORTANT: Run this script from the Codes folder!

# To run with default settings (uses your project structure):
result <- process_all_data()

# To create and use custom configuration:
# create_example_config("my_config.yaml")  # Edit this file as needed
# result <- process_all_data("my_config.yaml")

# Alternative: If running from different location, specify directories explicitly:
# config_with_paths <- list(
#   base_dir = "path/to/your/project/root",
#   split_data_dir = "path/to/your/project/root/Analysis/2. Split_Data", 
#   metadata_dir = "path/to/your/project/root/Data/Orthomosaics+shapefile/Metadata",
#   output_dir = "path/to/your/project/root/Analysis/4. Merged_data"
# )
# result <- process_all_data(config_with_paths)

# Print results
# if (result$success) {
#   cat("All files processed successfully!\n")
# } else {
#   cat("Some files failed to process. Check the log for details.\n")
# }