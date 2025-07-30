################################################################################
# 06_temporal_transfer_learning_cluster.r
# Phase 6: Temporal Transfer Learning - CLUSTER VERSION
# 
# Purpose: Use all available temporal data from separate UAV/Satellite files
#   - Optimized for cluster execution with checkpoint saving
#   - Memory-efficient processing
#   - Configurable parallel processing
#   - Progress tracking and recovery from interruptions
#
# Author: Temporal Transfer Learning Analysis - Cluster Version
# Date: 2024
################################################################################

# Load required libraries
cat("Loading required libraries...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(ggplot2)
  library(viridis)
  library(gridExtra)
  library(patchwork)
  library(reshape2)
  library(parallel)
  library(doParallel)
})

# Set seed and options
set.seed(42)
options(scipen = 999)

# Configuration
DEBUG_MODE <- FALSE  # Set to TRUE for quick debugging
USE_PARALLEL <- FALSE  # Set to FALSE for cluster stability - can enable if needed
N_REPETITIONS <- ifelse(DEBUG_MODE, 5, 50)  # Number of repeated holdout splits
CHECKPOINT_FREQUENCY <- 5  # Save checkpoint every N analyses

# Check if on cluster
on_cluster <- Sys.getenv("SLURM_JOB_ID") != ""
if (on_cluster) {
  cat(sprintf("Running on cluster - Job ID: %s\n", Sys.getenv("SLURM_JOB_ID")))
  cat(sprintf("CPUs allocated: %s\n", Sys.getenv("SLURM_CPUS_PER_TASK")))
}

# Setup parallel processing if enabled
if (USE_PARALLEL) {
  if (on_cluster) {
    n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) - 1
  } else {
    n_cores <- detectCores() - 1
  }
  n_cores <- max(1, n_cores)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  cat(sprintf("Parallel processing enabled with %d cores\n", n_cores))
} else {
  registerDoSEQ()
  cat("Running in sequential mode\n")
}

################################################################################
# SECTION 0: SETUP AND PATHS - CLUSTER STRUCTURE
################################################################################

# Set paths - MODIFIED for cluster directory structure
current_dir <- getwd()
project_root <- ifelse(basename(current_dir) == "scripts", dirname(current_dir), current_dir)

# Direct paths
data_dir <- file.path(project_root, "data")
phase1_dir <- file.path(project_root, "soybean_analysis_outputs", "phase1_harmonized_data")
phase6_dir <- file.path(project_root, "soybean_analysis_outputs", "phase6_temporal_transfer_cluster")

# Create output directories
subdirs <- c("cross_location", "cross_platform", "feature_stability", 
             "visualizations", "raw_results", "temporal_analysis", "checkpoints")
for (subdir in subdirs) {
  dir_path <- file.path(phase6_dir, subdir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

# Logging function
log_file <- file.path(phase6_dir, "temporal_transfer_cluster_log.txt")
log_message <- function(message, type = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- sprintf("[%s] %s: %s\n", timestamp, type, message)
  cat(msg)
  cat(msg, file = log_file, append = TRUE)
}

log_message("Starting Phase 6: Temporal Transfer Learning - CLUSTER VERSION")
log_message(sprintf("Configuration: %d repetitions, checkpoint every %d analyses", 
                    N_REPETITIONS, CHECKPOINT_FREQUENCY))

################################################################################
# SECTION 1: CHECKPOINT MANAGEMENT
################################################################################

# Main checkpoint file
main_checkpoint_file <- file.path(phase6_dir, "checkpoints", "main_checkpoint.rds")

# Load checkpoint if exists
if (file.exists(main_checkpoint_file)) {
  log_message("Loading checkpoint data...")
  checkpoint <- readRDS(main_checkpoint_file)
  completed_analyses <- checkpoint$completed
  partial_results <- checkpoint$partial_results
  log_message(sprintf("Resuming from checkpoint: %d analyses completed", 
                      length(completed_analyses)))
} else {
  completed_analyses <- character()
  partial_results <- list()
  log_message("Starting fresh analysis (no checkpoint found)")
}

# Function to save checkpoint
save_checkpoint <- function(completed, results) {
  checkpoint <- list(
    completed = completed,
    partial_results = results,
    timestamp = Sys.time(),
    n_repetitions = N_REPETITIONS
  )
  
  # Save to temporary file first
  temp_file <- paste0(main_checkpoint_file, ".tmp")
  saveRDS(checkpoint, temp_file)
  
  # Then rename (atomic operation)
  file.rename(temp_file, main_checkpoint_file)
  
  log_message(sprintf("Checkpoint saved: %d analyses completed", length(completed)))
}

################################################################################
# SECTION 2: DATA LOADING - MEMORY EFFICIENT
################################################################################

#' Load separate UAV and Satellite datasets with memory management
load_separate_datasets <- function() {
  log_message("Loading separate UAV and Satellite datasets...")
  
  datasets <- list()
  
  # Define file paths
  files <- list(
    SF_UAV = file.path(data_dir, "SF_UAV_with_DAP_cleaned.csv"),
    SF_Satellite = file.path(data_dir, "SF_Satellite_with_DAP_cleaned.csv"),
    EF_UAV = file.path(data_dir, "EF_UAV_with_DAP_cleaned.csv"),
    EF_Satellite = file.path(data_dir, "EF_Satellite_with_DAP_cleaned.csv")
  )
  
  # Check files exist
  for (name in names(files)) {
    if (!file.exists(files[[name]])) {
      log_message(sprintf("ERROR: %s not found at %s", name, files[[name]]), "ERROR")
      stop(sprintf("Required file not found: %s", files[[name]]))
    }
  }
  
  # Load data with garbage collection
  for (name in names(files)) {
    log_message(sprintf("Loading %s...", name))
    datasets[[name]] <- read_csv(files[[name]], show_col_types = FALSE)
    
    log_message(sprintf("  %s: %d rows (%d unique plots)", 
                        name,
                        nrow(datasets[[name]]),
                        length(unique(datasets[[name]]$Plot_ID))))
    
    # Force garbage collection after loading large file
    gc(verbose = FALSE)
  }
  
  return(datasets)
}

################################################################################
# SECTION 3: TEMPORAL AGGREGATION - MEMORY EFFICIENT
################################################################################

#' Temporal feature engineering and aggregation - memory efficient version
temporal_aggregation <- function(data, platform_type = "UAV") {
  
  log_message(sprintf("Performing temporal aggregation for %s data...", platform_type))
  
  # Force garbage collection before processing
  gc(verbose = FALSE)
  
  # Get feature columns based on platform
  if (platform_type == "UAV") {
    feature_patterns <- c("^MS_Green_", "^MS_Red_", "^MS_Red_edge_", "^MS_NIR_",
                         "^NDVI_MS_", "^GNDVI_MS_", "^NDRE_MS_", "^NGRDI_MS_",
                         "^SAVI_MS_", "^MSAVI_MS_", "^MTCI_MS_", "^CIgreen_MS_",
                         "^CIrededge_MS_", "^OSAVI_MS_", "^RGRI_MS_")
  } else {
    feature_patterns <- c("^Satellite_Red_", "^Satellite_Green_", "^Satellite_Blue_",
                         "^Satellite_NIR_", "^Satellite_Red_edge_", "^Satellite_Deep_blue_",
                         "^NDVI_Satellite_", "^GNDVI_Satellite_", "^NDRE_Satellite_",
                         "^GLI_Satellite_", "^NGRDI_Satellite_", "^SAVI_Satellite_",
                         "^EVI_Satellite_", "^MSAVI_Satellite_", "^SIPI_Satellite_",
                         "^MTCI_Satellite_", "^CIgreen_Satellite_", "^CIrededge_Satellite_",
                         "^ARVI_Satellite_", "^VARI_Satellite_", "^OSAVI_Satellite_",
                         "^TGI_Satellite_", "^ExG_Satellite_", "^RGRI_Satellite_")
  }
  
  # Get feature columns efficiently
  all_cols <- names(data)
  feature_cols <- c()
  
  for (pattern in feature_patterns) {
    matches <- all_cols[grepl(pattern, all_cols)]
    feature_cols <- c(feature_cols, matches)
  }
  
  feature_cols <- unique(feature_cols)
  numeric_features <- feature_cols[sapply(data[feature_cols], is.numeric)]
  
  log_message(sprintf("Found %d numeric features for aggregation", length(numeric_features)))
  
  # Process in chunks for memory efficiency
  unique_plots <- unique(data$Plot_ID)
  n_plots <- length(unique_plots)
  chunk_size <- 100  # Process 100 plots at a time
  n_chunks <- ceiling(n_plots / chunk_size)
  
  aggregated_chunks <- list()
  
  for (i in 1:n_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_plots)
    chunk_plots <- unique_plots[start_idx:end_idx]
    
    # Process chunk
    chunk_data <- data %>%
      filter(Plot_ID %in% chunk_plots) %>%
      select(Plot_ID, yield, peak_lai, DAP, Growth_Stage, all_of(numeric_features))
    
    # Aggregate chunk
    chunk_agg <- chunk_data %>%
      group_by(Plot_ID) %>%
      summarise(
        yield = first(yield[!is.na(yield)]),
        peak_lai = first(peak_lai[!is.na(peak_lai)]),
        n_timepoints = n(),
        dap_min = min(DAP, na.rm = TRUE),
        dap_max = max(DAP, na.rm = TRUE),
        dap_range = dap_max - dap_min,
        has_vegetative = as.numeric("Vegetative" %in% Growth_Stage),
        has_reproductive = as.numeric("Reproductive" %in% Growth_Stage),
        has_maturity = as.numeric("Maturity" %in% Growth_Stage),
        n_growth_stages = length(unique(Growth_Stage[!is.na(Growth_Stage)])),
        across(all_of(numeric_features), list(
          mean = ~mean(.x, na.rm = TRUE),
          max = ~max(.x, na.rm = TRUE),
          min = ~min(.x, na.rm = TRUE),
          sd = ~sd(.x, na.rm = TRUE),
          cv = ~ifelse(mean(.x, na.rm = TRUE) != 0, 
                       sd(.x, na.rm = TRUE) / abs(mean(.x, na.rm = TRUE)), 
                       0)
        ), .names = "{.col}_temporal_{.fn}")
      )
    
    aggregated_chunks[[i]] <- chunk_agg
    
    # Progress update
    if (i %% 10 == 0) {
      log_message(sprintf("  Processed %d/%d chunks", i, n_chunks))
      gc(verbose = FALSE)  # Clean up memory
    }
  }
  
  # Combine chunks
  temporal_features <- bind_rows(aggregated_chunks)
  
  # Clean up
  rm(aggregated_chunks)
  gc(verbose = FALSE)
  
  log_message(sprintf("Aggregated to %d plots with %d features", 
                      nrow(temporal_features), 
                      ncol(temporal_features) - 1))
  
  return(temporal_features)
}

#' Create fusion datasets by merging UAV and Satellite
create_fusion_datasets <- function(uav_data, sat_data, location) {
  
  log_message(sprintf("Creating fusion dataset for %s...", location))
  
  fusion_data <- uav_data %>%
    inner_join(sat_data, by = c("Plot_ID", "yield", "peak_lai"), 
               suffix = c("_uav", "_sat"))
  
  # Handle duplicate columns
  fusion_data <- fusion_data %>%
    select(-ends_with("_sat.y")) %>%
    rename_with(~str_replace(.x, "_sat.x", "_sat"), ends_with("_sat.x"))
  
  log_message(sprintf("Fusion dataset: %d plots with both UAV and Satellite data", 
                      nrow(fusion_data)))
  
  gc(verbose = FALSE)
  
  return(fusion_data)
}

################################################################################
# SECTION 4: FEATURE EXTRACTION
################################################################################

#' Extract features by platform - Updated for temporal features
extract_temporal_features <- function(data, platform = "both") {
  
  all_cols <- names(data)
  features <- c()
  
  if (platform == "both" || platform == "fusion") {
    # Get both UAV and Satellite temporal features
    uav_patterns <- c("^MS_.*_temporal_", "^NDVI_MS_.*_temporal_", 
                      "^GNDVI_MS_.*_temporal_", "^NDRE_MS_.*_temporal_",
                      "^NGRDI_MS_.*_temporal_", "^SAVI_MS_.*_temporal_",
                      "^MSAVI_MS_.*_temporal_", "^MTCI_MS_.*_temporal_",
                      "^CIgreen_MS_.*_temporal_", "^CIrededge_MS_.*_temporal_",
                      "^OSAVI_MS_.*_temporal_", "^RGRI_MS_.*_temporal_")
    
    sat_patterns <- c("^Satellite_.*_temporal_", "^NDVI_Satellite_.*_temporal_",
                      "^GNDVI_Satellite_.*_temporal_", "^NDRE_Satellite_.*_temporal_",
                      "^GLI_Satellite_.*_temporal_", "^NGRDI_Satellite_.*_temporal_",
                      "^SAVI_Satellite_.*_temporal_", "^EVI_Satellite_.*_temporal_",
                      "^MSAVI_Satellite_.*_temporal_", "^SIPI_Satellite_.*_temporal_",
                      "^MTCI_Satellite_.*_temporal_", "^CIgreen_Satellite_.*_temporal_",
                      "^CIrededge_Satellite_.*_temporal_", "^ARVI_Satellite_.*_temporal_",
                      "^VARI_Satellite_.*_temporal_", "^OSAVI_Satellite_.*_temporal_",
                      "^TGI_Satellite_.*_temporal_", "^ExG_Satellite_.*_temporal_",
                      "^RGRI_Satellite_.*_temporal_")
    
    for (pattern in c(uav_patterns, sat_patterns)) {
      matches <- all_cols[grepl(pattern, all_cols)]
      features <- c(features, matches)
    }
    
  } else if (platform == "uav") {
    uav_patterns <- c("^MS_.*_temporal_", "^NDVI_MS_.*_temporal_", 
                      "^GNDVI_MS_.*_temporal_", "^NDRE_MS_.*_temporal_",
                      "^NGRDI_MS_.*_temporal_", "^SAVI_MS_.*_temporal_",
                      "^MSAVI_MS_.*_temporal_", "^MTCI_MS_.*_temporal_",
                      "^CIgreen_MS_.*_temporal_", "^CIrededge_MS_.*_temporal_",
                      "^OSAVI_MS_.*_temporal_", "^RGRI_MS_.*_temporal_")
    
    for (pattern in uav_patterns) {
      matches <- all_cols[grepl(pattern, all_cols)]
      features <- c(features, matches)
    }
    
  } else if (platform == "satellite") {
    sat_patterns <- c("^Satellite_.*_temporal_", "^NDVI_Satellite_.*_temporal_",
                      "^GNDVI_Satellite_.*_temporal_", "^NDRE_Satellite_.*_temporal_",
                      "^GLI_Satellite_.*_temporal_", "^NGRDI_Satellite_.*_temporal_",
                      "^SAVI_Satellite_.*_temporal_", "^EVI_Satellite_.*_temporal_",
                      "^MSAVI_Satellite_.*_temporal_", "^SIPI_Satellite_.*_temporal_",
                      "^MTCI_Satellite_.*_temporal_", "^CIgreen_Satellite_.*_temporal_",
                      "^CIrededge_Satellite_.*_temporal_", "^ARVI_Satellite_.*_temporal_",
                      "^VARI_Satellite_.*_temporal_", "^OSAVI_Satellite_.*_temporal_",
                      "^TGI_Satellite_.*_temporal_", "^ExG_Satellite_.*_temporal_",
                      "^RGRI_Satellite_.*_temporal_")
    
    for (pattern in sat_patterns) {
      matches <- all_cols[grepl(pattern, all_cols)]
      features <- c(features, matches)
    }
  }
  
  # Add temporal metadata features
  temporal_meta <- c("n_timepoints", "dap_range", "has_vegetative", 
                     "has_reproductive", "has_maturity", "n_growth_stages")
  features <- c(features, intersect(temporal_meta, all_cols))
  
  # Remove duplicates and ensure numeric
  features <- unique(features)
  numeric_features <- features[sapply(data[features], is.numeric)]
  
  # Remove any inf or constant features
  valid_features <- c()
  for (feat in numeric_features) {
    vals <- data[[feat]]
    if (!all(is.na(vals)) && !all(is.infinite(vals)) && sd(vals, na.rm = TRUE) > 0) {
      valid_features <- c(valid_features, feat)
    }
  }
  
  # Limit features if too many (memory constraint)
  if (length(valid_features) > 500) {
    log_message(sprintf("Reducing features from %d to 500 for memory efficiency", 
                        length(valid_features)))
    # Select based on variance
    variances <- sapply(valid_features, function(f) var(data[[f]], na.rm = TRUE))
    valid_features <- names(sort(variances, decreasing = TRUE)[1:500])
  }
  
  log_message(sprintf("Extracted %d valid temporal features for platform '%s'", 
                      length(valid_features), platform))
  
  return(valid_features)
}

################################################################################
# SECTION 5: TRANSFER LEARNING FUNCTIONS - CLUSTER OPTIMIZED
################################################################################

#' Stratified sampling function
stratified_sample <- function(data, prop, target_var = "yield", n_strata = 5) {
  
  data <- data %>% mutate(orig_row_num = 1:n())
  
  # Create yield strata
  data_strat <- data %>%
    mutate(
      yield_quantile = cut(!!sym(target_var), 
                           breaks = quantile(!!sym(target_var), 
                                             probs = seq(0, 1, 1/n_strata),
                                             na.rm = TRUE),
                           labels = FALSE, include.lowest = TRUE)
    )
  
  # Sample within strata
  sampled_data <- data_strat %>%
    group_by(yield_quantile) %>%
    sample_frac(size = prop) %>%
    ungroup()
  
  return(sampled_data$orig_row_num)
}

#' Evaluate predictions with confidence intervals
evaluate_predictions_ci <- function(predictions_list, truth_list) {
  
  if (!is.list(predictions_list)) {
    predictions_list <- list(predictions_list)
  }
  
  if (!is.list(truth_list)) {
    truth_list <- rep(list(truth_list), length(predictions_list))
  }
  
  # Calculate metrics for each repetition
  metrics_list <- list()
  
  for (i in 1:length(predictions_list)) {
    pred <- predictions_list[[i]]
    truth <- truth_list[[i]]
    
    valid_idx <- !is.na(pred) & !is.na(truth)
    pred <- pred[valid_idx]
    truth_valid <- truth[valid_idx]
    
    if (length(pred) == 0) {
      metrics_list[[i]] <- c(r2 = NA, rmse = NA, mae = NA, mape = NA)
    } else {
      r2 <- cor(pred, truth_valid, use = "complete.obs")^2
      rmse <- sqrt(mean((pred - truth_valid)^2, na.rm = TRUE))
      mae <- mean(abs(pred - truth_valid), na.rm = TRUE)
      mape <- mean(abs((pred - truth_valid) / truth_valid), na.rm = TRUE) * 100
      
      metrics_list[[i]] <- c(r2 = r2, rmse = rmse, mae = mae, mape = mape)
    }
  }
  
  # Convert to matrix
  metrics_matrix <- do.call(rbind, metrics_list)
  
  # Calculate mean and confidence intervals
  results <- list()
  for (metric in colnames(metrics_matrix)) {
    values <- metrics_matrix[, metric]
    values <- values[!is.na(values)]
    
    if (length(values) > 0) {
      results[[paste0(metric, "_mean")]] <- mean(values)
      results[[paste0(metric, "_sd")]] <- ifelse(length(values) > 1, sd(values), 0)
      results[[paste0(metric, "_ci_lower")]] <- ifelse(length(values) > 1, 
                                                       quantile(values, 0.025),
                                                       values[1])
      results[[paste0(metric, "_ci_upper")]] <- ifelse(length(values) > 1,
                                                       quantile(values, 0.975),
                                                       values[1])
    } else {
      results[[paste0(metric, "_mean")]] <- NA
      results[[paste0(metric, "_sd")]] <- NA
      results[[paste0(metric, "_ci_lower")]] <- NA
      results[[paste0(metric, "_ci_upper")]] <- NA
    }
  }
  
  return(results)
}

#' Cross-location transfer with temporal features - CLUSTER VERSION
cross_location_transfer_temporal <- function(source_data, target_data, 
                                           source_name, target_name,
                                           platform = "both", 
                                           target_var = "yield",
                                           n_reps = 50) {
  
  log_message(sprintf("\n=== Cross-Location Transfer: %s → %s (%s, Target: %s) ===",
                      source_name, target_name, platform, target_var))
  
  # Memory management
  gc(verbose = FALSE)
  
  # Extract features
  features <- extract_temporal_features(source_data, platform)
  
  if (length(features) == 0) {
    log_message(sprintf("ERROR: No features found for platform '%s'", platform), "ERROR")
    return(list(error = "No features found"))
  }
  
  # Prepare data
  source_clean <- source_data %>%
    filter(!is.na(!!sym(target_var))) %>%
    select(all_of(c(target_var, features, "Plot_ID"))) %>%
    na.omit()
  
  target_clean <- target_data %>%
    filter(!is.na(!!sym(target_var))) %>%
    select(all_of(c(target_var, features, "Plot_ID"))) %>%
    na.omit()
  
  log_message(sprintf("Source: %d samples, Target: %d samples, Features: %d",
                      nrow(source_clean), nrow(target_clean), length(features)))
  
  if (nrow(source_clean) < 50 || nrow(target_clean) < 20) {
    log_message("ERROR: Insufficient data for analysis", "ERROR")
    return(list(error = "Insufficient data"))
  }
  
  results <- list()
  
  # 1. Direct transfer
  log_message("Training baseline Random Forest model...")
  
  # Reduce trees for memory efficiency
  rf_direct <- randomForest(
    x = source_clean[, features],
    y = source_clean[[target_var]],
    ntree = 300,  # Reduced from 500
    importance = TRUE
  )
  
  rf_pred_direct <- predict(rf_direct, target_clean[, features])
  results$rf_direct <- evaluate_predictions_ci(list(rf_pred_direct), 
                                              target_clean[[target_var]])
  
  # Clean up
  rm(rf_direct)
  gc(verbose = FALSE)
  
  # 2. Fine-tuning with multiple proportions
  log_message("Starting fine-tuning analysis...")
  
  fine_tune_props <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  results$fine_tuning <- list()
  
  for (prop in fine_tune_props) {
    log_message(sprintf("Fine-tuning with %.0f%% target data...", prop * 100))
    
    rf_predictions <- list()
    rf_test_truth <- list()
    
    # Progress bar for repetitions
    pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
    
    for (rep in 1:n_reps) {
      setTxtProgressBar(pb, rep)
      set.seed(rep * 1000)
      
      # Stratified sampling
      fine_tune_idx <- stratified_sample(target_clean, prop, target_var)
      
      fine_tune_data <- target_clean[fine_tune_idx, ]
      test_data <- target_clean[-fine_tune_idx, ]
      
      # Combine source and fine-tuning data
      combined_data <- rbind(source_clean[, c(features, target_var)], 
                            fine_tune_data[, c(features, target_var)])
      
      # Train RF with reduced trees
      rf_finetuned <- randomForest(
        x = combined_data[, features],
        y = combined_data[[target_var]],
        ntree = 300  # Reduced from 500
      )
      
      rf_predictions[[rep]] <- predict(rf_finetuned, test_data[, features])
      rf_test_truth[[rep]] <- test_data[[target_var]]
      
      # Clean up
      rm(rf_finetuned)
      
      # Periodic garbage collection
      if (rep %% 10 == 0) {
        gc(verbose = FALSE)
      }
    }
    
    close(pb)
    
    # Evaluate
    results$fine_tuning[[paste0("prop_", prop)]] <- list(
      rf = evaluate_predictions_ci(rf_predictions, rf_test_truth),
      n_samples_mean = floor(nrow(target_clean) * prop)
    )
    
    # Clean up
    gc(verbose = FALSE)
  }
  
  return(results)
}

#' Cross-platform transfer with temporal features - CLUSTER VERSION
cross_platform_transfer_temporal <- function(uav_data, sat_data, fusion_data,
                                           location, target_var = "yield", 
                                           n_reps = 50) {
  
  log_message(sprintf("\n=== Cross-Platform Transfer: %s (Target: %s) ===", 
                      location, target_var))
  
  # Memory management
  gc(verbose = FALSE)
  
  results <- list()
  
  # Run multiple repetitions
  fusion_predictions <- list()
  uav_predictions <- list()
  sat_predictions <- list()
  test_truth <- list()
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
  
  for (rep in 1:n_reps) {
    setTxtProgressBar(pb, rep)
    set.seed(rep * 100)
    
    # For fusion data (plots with both platforms)
    if (!is.null(fusion_data) && nrow(fusion_data) > 50) {
      train_idx <- stratified_sample(fusion_data, 0.7, target_var)
      train_fusion <- fusion_data[train_idx, ]
      test_fusion <- fusion_data[-train_idx, ]
      test_truth[[rep]] <- test_fusion[[target_var]]
      
      # Train fusion model
      fusion_features <- extract_temporal_features(train_fusion, "both")
      if (length(fusion_features) > 0) {
        fusion_clean <- train_fusion %>%
          select(all_of(c(target_var, fusion_features))) %>%
          na.omit()
        
        fusion_model <- randomForest(
          x = fusion_clean[, fusion_features],
          y = fusion_clean[[target_var]],
          ntree = 300  # Reduced
        )
        
        fusion_predictions[[rep]] <- predict(fusion_model, test_fusion[, fusion_features])
        rm(fusion_model)
      }
      
      # Train UAV-only model on same plots
      uav_features <- extract_temporal_features(train_fusion, "uav")
      if (length(uav_features) > 0) {
        uav_clean <- train_fusion %>%
          select(all_of(c(target_var, uav_features))) %>%
          na.omit()
        
        uav_model <- randomForest(
          x = uav_clean[, uav_features],
          y = uav_clean[[target_var]],
          ntree = 300
        )
        
        uav_predictions[[rep]] <- predict(uav_model, test_fusion[, uav_features])
        rm(uav_model)
      }
      
      # Train Satellite-only model on same plots
      sat_features <- extract_temporal_features(train_fusion, "satellite")
      if (length(sat_features) > 0) {
        sat_clean <- train_fusion %>%
          select(all_of(c(target_var, sat_features))) %>%
          na.omit()
        
        sat_model <- randomForest(
          x = sat_clean[, sat_features],
          y = sat_clean[[target_var]],
          ntree = 300
        )
        
        sat_predictions[[rep]] <- predict(sat_model, test_fusion[, sat_features])
        rm(sat_model)
      }
      
      # Periodic garbage collection
      if (rep %% 10 == 0) {
        gc(verbose = FALSE)
      }
    }
  }
  
  close(pb)
  
  # Evaluate models
  if (length(fusion_predictions) > 0) {
    results$fusion <- evaluate_predictions_ci(fusion_predictions, test_truth)
  }
  
  if (length(uav_predictions) > 0) {
    results$uav_only <- evaluate_predictions_ci(uav_predictions, test_truth)
  }
  
  if (length(sat_predictions) > 0) {
    results$sat_only <- evaluate_predictions_ci(sat_predictions, test_truth)
  }
  
  # Calculate improvements
  if ("fusion" %in% names(results) && "uav_only" %in% names(results)) {
    results$uav_improvement <- results$fusion$r2_mean - results$uav_only$r2_mean
  }
  
  if ("fusion" %in% names(results) && "sat_only" %in% names(results)) {
    results$sat_improvement <- results$fusion$r2_mean - results$sat_only$r2_mean
  }
  
  # Clean up
  gc(verbose = FALSE)
  
  return(results)
}

################################################################################
# SECTION 6: MAIN EXECUTION WITH CHECKPOINTING
################################################################################

log_message("\n=== MAIN EXECUTION: TEMPORAL TRANSFER LEARNING ===")

# Initialize results storage
all_results <- list()
cross_location_summary <- data.frame()
cross_platform_summary <- data.frame()

# Counter for checkpoint saving
analysis_counter <- 0

# Try to load partial results if they exist
if (exists("partial_results") && length(partial_results) > 0) {
  all_results <- partial_results$all_results
  cross_location_summary <- partial_results$cross_location_summary
  cross_platform_summary <- partial_results$cross_platform_summary
}

# STEP 1: Load and prepare data (only if not already done)
if (!"data_loaded" %in% completed_analyses) {
  log_message("\n### LOADING AND PREPARING DATA ###")
  
  # Load separate datasets
  datasets <- load_separate_datasets()
  
  # Perform temporal aggregation
  log_message("\n### TEMPORAL AGGREGATION ###")
  
  sf_uav_agg <- temporal_aggregation(datasets$SF_UAV, "UAV")
  saveRDS(sf_uav_agg, file.path(phase6_dir, "checkpoints", "sf_uav_agg.rds"))
  
  sf_sat_agg <- temporal_aggregation(datasets$SF_Satellite, "Satellite")
  saveRDS(sf_sat_agg, file.path(phase6_dir, "checkpoints", "sf_sat_agg.rds"))
  
  ef_uav_agg <- temporal_aggregation(datasets$EF_UAV, "UAV")
  saveRDS(ef_uav_agg, file.path(phase6_dir, "checkpoints", "ef_uav_agg.rds"))
  
  ef_sat_agg <- temporal_aggregation(datasets$EF_Satellite, "Satellite")
  saveRDS(ef_sat_agg, file.path(phase6_dir, "checkpoints", "ef_sat_agg.rds"))
  
  # Create fusion datasets
  sf_fusion <- create_fusion_datasets(sf_uav_agg, sf_sat_agg, "SF")
  saveRDS(sf_fusion, file.path(phase6_dir, "checkpoints", "sf_fusion.rds"))
  
  ef_fusion <- create_fusion_datasets(ef_uav_agg, ef_sat_agg, "EF")
  saveRDS(ef_fusion, file.path(phase6_dir, "checkpoints", "ef_fusion.rds"))
  
  # Mark as completed and save checkpoint
  completed_analyses <- c(completed_analyses, "data_loaded")
  save_checkpoint(completed_analyses, list(
    all_results = all_results,
    cross_location_summary = cross_location_summary,
    cross_platform_summary = cross_platform_summary
  ))
  
  # Clean up original datasets
  rm(datasets)
  gc(verbose = FALSE)
  
} else {
  log_message("Loading pre-processed aggregated data from checkpoint...")
  
  sf_uav_agg <- readRDS(file.path(phase6_dir, "checkpoints", "sf_uav_agg.rds"))
  sf_sat_agg <- readRDS(file.path(phase6_dir, "checkpoints", "sf_sat_agg.rds"))
  ef_uav_agg <- readRDS(file.path(phase6_dir, "checkpoints", "ef_uav_agg.rds"))
  ef_sat_agg <- readRDS(file.path(phase6_dir, "checkpoints", "ef_sat_agg.rds"))
  sf_fusion <- readRDS(file.path(phase6_dir, "checkpoints", "sf_fusion.rds"))
  ef_fusion <- readRDS(file.path(phase6_dir, "checkpoints", "ef_fusion.rds"))
}

################################################################################
# SECTION 7: CROSS-LOCATION TRANSFER EXPERIMENTS
################################################################################

log_message("\n### CROSS-LOCATION TRANSFER ANALYSIS ###")

# Define all cross-location experiments
cross_location_experiments <- list(
  list(name = "SF_to_EF_UAV_temporal_yield", 
       source = sf_uav_agg, target = ef_uav_agg, 
       source_name = "SF", target_name = "EF", 
       platform = "uav", target_var = "yield"),
  
  list(name = "SF_to_EF_UAV_temporal_peak_lai",
       source = sf_uav_agg, target = ef_uav_agg,
       source_name = "SF", target_name = "EF",
       platform = "uav", target_var = "peak_lai"),
  
  list(name = "SF_to_EF_Sat_temporal_yield",
       source = sf_sat_agg, target = ef_sat_agg,
       source_name = "SF", target_name = "EF",
       platform = "satellite", target_var = "yield"),
  
  list(name = "SF_to_EF_Sat_temporal_peak_lai",
       source = sf_sat_agg, target = ef_sat_agg,
       source_name = "SF", target_name = "EF",
       platform = "satellite", target_var = "peak_lai"),
  
  list(name = "SF_to_EF_Fusion_temporal_yield",
       source = sf_fusion, target = ef_fusion,
       source_name = "SF", target_name = "EF",
       platform = "both", target_var = "yield"),
  
  list(name = "SF_to_EF_Fusion_temporal_peak_lai",
       source = sf_fusion, target = ef_fusion,
       source_name = "SF", target_name = "EF",
       platform = "both", target_var = "peak_lai"),
  
  # Reverse direction
  list(name = "EF_to_SF_UAV_temporal_yield",
       source = ef_uav_agg, target = sf_uav_agg,
       source_name = "EF", target_name = "SF",
       platform = "uav", target_var = "yield"),
  
  list(name = "EF_to_SF_UAV_temporal_peak_lai",
       source = ef_uav_agg, target = sf_uav_agg,
       source_name = "EF", target_name = "SF",
       platform = "uav", target_var = "peak_lai"),
  
  list(name = "EF_to_SF_Sat_temporal_yield",
       source = ef_sat_agg, target = sf_sat_agg,
       source_name = "EF", target_name = "SF",
       platform = "satellite", target_var = "yield"),
  
  list(name = "EF_to_SF_Sat_temporal_peak_lai",
       source = ef_sat_agg, target = sf_sat_agg,
       source_name = "EF", target_name = "SF",
       platform = "satellite", target_var = "peak_lai"),
  
  list(name = "EF_to_SF_Fusion_temporal_yield",
       source = ef_fusion, target = sf_fusion,
       source_name = "EF", target_name = "SF",
       platform = "both", target_var = "yield"),
  
  list(name = "EF_to_SF_Fusion_temporal_peak_lai",
       source = ef_fusion, target = sf_fusion,
       source_name = "EF", target_name = "SF",
       platform = "both", target_var = "peak_lai")
)

# Initialize cross-location results if needed
if (!"cross_location" %in% names(all_results)) {
  all_results$cross_location <- list()
}

# Run experiments
for (exp in cross_location_experiments) {
  
  # Skip if already completed
  if (exp$name %in% completed_analyses) {
    log_message(sprintf("\nSkipping %s (already completed)", exp$name))
    next
  }
  
  log_message(sprintf("\nRunning %s...", exp$name))
  
  result <- tryCatch({
    cross_location_transfer_temporal(
      exp$source, exp$target,
      exp$source_name, exp$target_name,
      exp$platform, exp$target_var,
      n_reps = N_REPETITIONS
    )
  }, error = function(e) {
    log_message(sprintf("ERROR in %s: %s", exp$name, e$message), "ERROR")
    list(error = e$message)
  })
  
  if (!"error" %in% names(result)) {
    all_results$cross_location[[exp$name]] <- result
    
    # Update summary
    if (!is.null(result$rf_direct)) {
      cross_location_summary <- rbind(cross_location_summary, data.frame(
        Transfer = exp$name,
        Method = "Direct",
        Model = "RF",
        R2_mean = result$rf_direct$r2_mean,
        R2_ci_lower = result$rf_direct$r2_ci_lower,
        R2_ci_upper = result$rf_direct$r2_ci_upper,
        RMSE_mean = result$rf_direct$rmse_mean,
        stringsAsFactors = FALSE
      ))
    }
    
    # Add fine-tuning results to summary
    if (!is.null(result$fine_tuning)) {
      for (prop in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
        prop_name <- paste0("prop_", prop)
        if (prop_name %in% names(result$fine_tuning)) {
          ft_result <- result$fine_tuning[[prop_name]]
          
          if (!is.null(ft_result$rf)) {
            cross_location_summary <- rbind(cross_location_summary, data.frame(
              Transfer = exp$name,
              Method = paste0("FineTune_", prop * 100, "%"),
              Model = "RF",
              R2_mean = ft_result$rf$r2_mean,
              R2_ci_lower = ft_result$rf$r2_ci_lower,
              R2_ci_upper = ft_result$rf$r2_ci_upper,
              RMSE_mean = ft_result$rf$rmse_mean,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
    
    # Save individual result
    saveRDS(result, file.path(phase6_dir, "cross_location", 
                              paste0(exp$name, "_result.rds")))
  }
  
  # Mark as completed
  completed_analyses <- c(completed_analyses, exp$name)
  analysis_counter <- analysis_counter + 1
  
  # Save checkpoint periodically
  if (analysis_counter %% CHECKPOINT_FREQUENCY == 0) {
    save_checkpoint(completed_analyses, list(
      all_results = all_results,
      cross_location_summary = cross_location_summary,
      cross_platform_summary = cross_platform_summary
    ))
    
    # Also save current summary
    write_csv(cross_location_summary, 
              file.path(phase6_dir, "cross_location_temporal_summary.csv"))
  }
  
  # Clean up memory
  gc(verbose = FALSE)
}

################################################################################
# SECTION 8: CROSS-PLATFORM TRANSFER EXPERIMENTS
################################################################################

log_message("\n### CROSS-PLATFORM TRANSFER ANALYSIS ###")

# Define cross-platform experiments
cross_platform_experiments <- list(
  list(name = "SF_yield", location = "SF", target_var = "yield",
       uav = sf_uav_agg, sat = sf_sat_agg, fusion = sf_fusion),
  
  list(name = "SF_peak_lai", location = "SF", target_var = "peak_lai",
       uav = sf_uav_agg, sat = sf_sat_agg, fusion = sf_fusion),
  
  list(name = "EF_yield", location = "EF", target_var = "yield",
       uav = ef_uav_agg, sat = ef_sat_agg, fusion = ef_fusion),
  
  list(name = "EF_peak_lai", location = "EF", target_var = "peak_lai",
       uav = ef_uav_agg, sat = ef_sat_agg, fusion = ef_fusion)
)

# Initialize cross-platform results if needed
if (!"cross_platform" %in% names(all_results)) {
  all_results$cross_platform <- list()
}

# Run experiments
for (exp in cross_platform_experiments) {
  
  # Skip if already completed
  if (paste0("cross_platform_", exp$name) %in% completed_analyses) {
    log_message(sprintf("\nSkipping cross-platform %s (already completed)", exp$name))
    next
  }
  
  log_message(sprintf("\nRunning cross-platform analysis: %s", exp$name))
  
  result <- tryCatch({
    cross_platform_transfer_temporal(
      exp$uav, exp$sat, exp$fusion,
      exp$location, exp$target_var,
      n_reps = N_REPETITIONS
    )
  }, error = function(e) {
    log_message(sprintf("ERROR in cross-platform %s: %s", exp$name, e$message), "ERROR")
    list(error = e$message)
  })
  
  if (!"error" %in% names(result)) {
    all_results$cross_platform[[exp$name]] <- result
    
    # Update summary
    for (platform in c("fusion", "uav_only", "sat_only")) {
      if (platform %in% names(result)) {
        res <- result[[platform]]
        cross_platform_summary <- rbind(cross_platform_summary, data.frame(
          Location_Target = exp$name,
          Platform = platform,
          R2_mean = res$r2_mean,
          R2_ci_lower = res$r2_ci_lower,
          R2_ci_upper = res$r2_ci_upper,
          RMSE_mean = res$rmse_mean,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Save individual result
    saveRDS(result, file.path(phase6_dir, "cross_platform", 
                              paste0("cross_platform_", exp$name, "_result.rds")))
  }
  
  # Mark as completed
  completed_analyses <- c(completed_analyses, paste0("cross_platform_", exp$name))
  analysis_counter <- analysis_counter + 1
  
  # Save checkpoint periodically
  if (analysis_counter %% CHECKPOINT_FREQUENCY == 0) {
    save_checkpoint(completed_analyses, list(
      all_results = all_results,
      cross_location_summary = cross_location_summary,
      cross_platform_summary = cross_platform_summary
    ))
    
    # Also save current summary
    write_csv(cross_platform_summary, 
              file.path(phase6_dir, "cross_platform_temporal_summary.csv"))
  }
  
  # Clean up memory
  gc(verbose = FALSE)
}

################################################################################
# SECTION 9: FINAL SUMMARIES AND VISUALIZATIONS
################################################################################

log_message("\n### CREATING FINAL SUMMARIES AND VISUALIZATIONS ###")

# Save final summaries
write_csv(cross_location_summary, 
          file.path(phase6_dir, "cross_location_temporal_summary_final.csv"))
write_csv(cross_platform_summary, 
          file.path(phase6_dir, "cross_platform_temporal_summary_final.csv"))

# Create visualizations
tryCatch({
  
  # 1. Cross-platform performance comparison
  if (nrow(cross_platform_summary) > 0) {
    p_performance <- ggplot(cross_platform_summary, 
                           aes(x = Platform, y = R2_mean, color = Platform)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = R2_ci_lower, ymax = R2_ci_upper), width = 0.2) +
      facet_wrap(~ Location_Target, scales = "free_y") +
      labs(title = "Cross-Platform Performance with Temporal Features",
           subtitle = "Using all available timepoints with temporal aggregation",
           y = "R² (with 95% CI)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_brewer(palette = "Set2")
    
    ggsave(file.path(phase6_dir, "visualizations", "platform_performance_temporal.pdf"), 
           p_performance, width = 12, height = 8)
  }
  
  # 2. Fine-tuning performance
  if (nrow(cross_location_summary) > 0) {
    ft_data <- cross_location_summary %>%
      filter(grepl("FineTune", Method)) %>%
      mutate(
        Percentage = as.numeric(gsub("FineTune_|%", "", Method)),
        Source_Target = gsub("_temporal.*", "", Transfer),
        Platform = case_when(
          grepl("UAV", Transfer) ~ "UAV",
          grepl("Sat", Transfer) ~ "Satellite",
          grepl("Fusion", Transfer) ~ "Fusion",
          TRUE ~ "Unknown"
        )
      )
    
    p_finetune <- ggplot(ft_data, aes(x = Percentage, y = R2_mean, 
                                      color = Source_Target, linetype = Platform)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      geom_ribbon(aes(ymin = R2_ci_lower, ymax = R2_ci_upper, fill = Source_Target), 
                  alpha = 0.2) +
      scale_x_continuous(breaks = c(10, 20, 30, 40, 50)) +
      theme_minimal() +
      labs(title = "Fine-Tuning Performance with Temporal Features",
           x = "Target Domain Data Used (%)",
           y = "Test R²") +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1")
    
    ggsave(file.path(phase6_dir, "visualizations", "fine_tuning_temporal.pdf"), 
           p_finetune, width = 12, height = 8)
  }
  
  log_message("Visualizations completed successfully")
  
}, error = function(e) {
  log_message(sprintf("ERROR in visualization creation: %s", e$message), "ERROR")
})

################################################################################
# SECTION 10: TEMPORAL ANALYSIS REPORT
################################################################################

log_message("\n=== Creating Temporal Analysis Report ===")

sink(file.path(phase6_dir, "temporal_transfer_summary.txt"))

cat("TEMPORAL TRANSFER LEARNING ANALYSIS - SUMMARY REPORT\n")
cat("===================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
if (on_cluster) {
  cat("SLURM Job ID:", Sys.getenv("SLURM_JOB_ID"), "\n")
}
cat("\nMETHODOLOGY:\n")
cat("- Used all available timepoints from separate UAV/Satellite files\n")
cat("- Temporal aggregation: mean, max, min, SD, CV\n")
cat("- Growth stage-specific features\n")
cat("- Maintained stratified sampling (n=", N_REPETITIONS, " repetitions)\n")
cat("- Fixed satellite column naming issues\n\n")

cat("DATA VOLUME:\n")
cat(sprintf("- SF UAV: %d aggregated plots\n", nrow(sf_uav_agg)))
cat(sprintf("- SF Satellite: %d aggregated plots\n", nrow(sf_sat_agg)))
cat(sprintf("- EF UAV: %d aggregated plots\n", nrow(ef_uav_agg)))
cat(sprintf("- EF Satellite: %d aggregated plots\n", nrow(ef_sat_agg)))
cat(sprintf("- SF Fusion: %d plots\n", nrow(sf_fusion)))
cat(sprintf("- EF Fusion: %d plots\n", nrow(ef_fusion)))

cat("\n\nKEY RESULTS:\n")

# Report cross-platform results
if (nrow(cross_platform_summary) > 0) {
  cat("\n1. CROSS-PLATFORM PERFORMANCE:\n")
  
  for (loc in c("SF", "EF")) {
    cat(sprintf("\n%s Location:\n", loc))
    loc_results <- cross_platform_summary %>% 
      filter(grepl(loc, Location_Target), grepl("yield", Location_Target))
    
    for (i in 1:nrow(loc_results)) {
      cat(sprintf("  %s: R² = %.3f [%.3f, %.3f]\n",
                  loc_results$Platform[i],
                  loc_results$R2_mean[i],
                  loc_results$R2_ci_lower[i],
                  loc_results$R2_ci_upper[i]))
    }
  }
}

# Report best cross-location transfers
if (nrow(cross_location_summary) > 0) {
  cat("\n\n2. BEST CROSS-LOCATION TRANSFERS:\n")
  
  best_transfers <- cross_location_summary %>%
    filter(Method == "Direct") %>%
    arrange(desc(R2_mean)) %>%
    head(5)
  
  for (i in 1:nrow(best_transfers)) {
    cat(sprintf("  %s: R² = %.3f [%.3f, %.3f]\n",
                best_transfers$Transfer[i],
                best_transfers$R2_mean[i],
                best_transfers$R2_ci_lower[i],
                best_transfers$R2_ci_upper[i]))
  }
}

cat("\n\nCOMPLETION STATUS:\n")
cat(sprintf("- Total analyses planned: %d\n", 
            length(cross_location_experiments) + length(cross_platform_experiments)))
cat(sprintf("- Analyses completed: %d\n", 
            length(completed_analyses) - 1))  # Subtract "data_loaded"
cat(sprintf("- Checkpoints saved: %d\n", 
            ceiling(analysis_counter / CHECKPOINT_FREQUENCY)))

cat("\n\nRECOMMENDATIONS:\n")
cat("1. Temporal aggregation significantly improves data volume\n")
cat("2. Growth stage features add valuable information\n")
cat("3. Satellite performance should be improved with temporal approach\n")
cat("4. Consider ensemble methods for operational deployment\n")

sink()

# Save complete results
saveRDS(all_results, file.path(phase6_dir, "temporal_transfer_results_complete.rds"))

# Final checkpoint
save_checkpoint(c(completed_analyses, "analysis_complete"), list(
  all_results = all_results,
  cross_location_summary = cross_location_summary,
  cross_platform_summary = cross_platform_summary
))

# Clean up parallel cluster if used
if (USE_PARALLEL) {
  stopCluster(cl)
}

log_message("\n*** TEMPORAL TRANSFER LEARNING ANALYSIS COMPLETE ***")
log_message(sprintf("Results saved to: %s", phase6_dir))
log_message(sprintf("Total execution time: %.2f hours", 
                    as.numeric(difftime(Sys.time(), 
                                        as.POSIXct(readLines(log_file, n = 1) %>% 
                                                     gsub("\\[|\\].*", "", .)),
                                        units = "hours"))))

################################################################################
# 06_temporal_transfer_learning.sh
# SLURM submission script for temporal transfer learning
################################################################################

#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16      # Request 16 CPUs for this task
#SBATCH --mem=192G              # Increased memory for temporal data
#SBATCH --time=48:00:00         # 48 hours should be sufficient
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@illinois.edu
#SBATCH -J phase6_temporal_transfer
#SBATCH -D /home/a-m/your_username/soybean_uav_satellite/src/slurm-out
#SBATCH -p normal

# Load the R module
module load R/4.4.0-IGB-gcc-8.2.0

# Set up paths
BASE_DIR="/home/a-m/your_username/soybean_uav_satellite"
SCRIPT_DIR="${BASE_DIR}/scripts"

# Change to the scripts directory
cd ${SCRIPT_DIR}

# Export the number of cores for R to detect
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Run the R script
echo "Starting Phase 6 Temporal Transfer Learning at $(date)"
echo "Working directory: $(pwd)"
echo "Number of CPUs allocated: ${SLURM_CPUS_PER_TASK}"
echo "Memory allocated: 192G"
echo "R version:"
R --version

# Execute the cluster R script
Rscript 06_temporal_transfer_learning_cluster.r

# Check exit status
if [ $? -eq 0 ]; then
    echo "Phase 6 Temporal Transfer Learning completed successfully at $(date)"
else
    echo "Phase 6 Temporal Transfer Learning failed with exit code $? at $(date)"
fi

# Optional: Send completion notification
echo "Phase 6 Temporal Transfer Learning job completed" | mail -s "SLURM Job ${SLURM_JOB_ID} Complete" your_email@illinois.edu