################################################################################
# Phase 5 Simplified: Multi-Platform Complementarity Analysis for Yield Prediction
#                     WITHOUT STACKED FUSION
# 
# Purpose: Simplified version that focuses on:
#   - Platform complementarity analysis
#   - Individual platform performance comparison
#   - Feature importance analysis
#   - Statistical testing
#   - NO stacked fusion components
#
# Author: Simplified Multi-Platform Analysis for Yield
# Date: 2024
################################################################################

# Load required libraries
cat("Loading required libraries...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(xgboost)
  library(corrplot)
  library(ggplot2)
  library(viridis)
  library(gridExtra)
  library(patchwork)
  library(reshape2)
  library(parallel)
  library(doParallel)
  library(plotly)
  library(boot)
})

# Set seed and options
set.seed(42)
options(scipen = 999)

# Configuration
DEBUG_MODE <- FALSE  # Set to TRUE for quick debugging
USE_PARALLEL <- TRUE  # Enable parallel processing
N_BOOTSTRAP <- ifelse(DEBUG_MODE, 10, 100)  # Bootstrap iterations
CONFIDENCE_LEVEL <- 0.95
SAVE_INTERVAL <- 10  # Save checkpoint every 10 iterations

# Check if on cluster
on_cluster <- Sys.getenv("SLURM_JOB_ID") != ""
if (on_cluster) {
  cat(sprintf("Running on cluster - Job ID: %s\n", Sys.getenv("SLURM_JOB_ID")))
  cat(sprintf("Allocated CPUs: %s\n", Sys.getenv("SLURM_CPUS_PER_TASK")))
}

# Start overall timer
overall_start_time <- Sys.time()

################################################################################
# SECTION 0: SETUP AND PATHS
################################################################################

# Set paths
current_dir <- getwd()
if (basename(current_dir) == "scripts") {
  project_root <- dirname(current_dir)
} else {
  project_root <- current_dir
}

# Direct paths
data_dir <- file.path(project_root, "data")
phase5_dir <- file.path(project_root, "soybean_analysis_outputs", "phase5_multiplatform_complementarity_yield")

# Create directories
subdirs <- c("complementarity_results", "feature_analysis", "model_performance", 
             "visualizations", "bootstrap_analysis", "confidence_intervals", 
             "summary_tables", "statistical_tests", "checkpoints", 
             "time_logs", "error_logs", "raw_outputs")
for (subdir in subdirs) {
  dir_path <- file.path(phase5_dir, subdir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

# Create location-specific directories
for (loc in c("SF", "EF")) {
  loc_subdirs <- c("platform_comparison", "feature_importance", "predictions", 
                   "bootstrap_results")
  for (subdir in loc_subdirs) {
    dir_path <- file.path(phase5_dir, "raw_outputs", loc, subdir)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
}

# Enhanced logging function
log_file <- file.path(phase5_dir, "complementarity_yield_analysis_log.txt")
error_log_file <- file.path(phase5_dir, "error_logs", "detailed_errors.txt")
time_log_file <- file.path(phase5_dir, "time_logs", "execution_times.csv")

log_message <- function(message, type = "INFO", log_time = FALSE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- sprintf("[%s] %s: %s\n", timestamp, type, message)
  cat(msg)
  cat(msg, file = log_file, append = TRUE)
  
  # Log execution time if requested
  if (log_time) {
    time_entry <- data.frame(
      timestamp = timestamp,
      message = message,
      elapsed_mins = as.numeric(difftime(Sys.time(), overall_start_time, units = "mins"))
    )
    if (file.exists(time_log_file)) {
      write.table(time_entry, time_log_file, append = TRUE, 
                  col.names = FALSE, row.names = FALSE, sep = ",")
    } else {
      write.csv(time_entry, time_log_file, row.names = FALSE)
    }
  }
}

log_message("Starting Simplified Multi-Platform Complementarity Analysis for YIELD")
log_message(sprintf("Configuration: %d bootstrap iterations, %.0f%% CI",
                   N_BOOTSTRAP, CONFIDENCE_LEVEL * 100))

# Setup parallel processing
if (USE_PARALLEL) {
  n_cores <- ifelse(on_cluster, 
                    min(8, as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", detectCores() - 1))),
                    min(4, detectCores() - 1))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  log_message(sprintf("Parallel processing enabled with %d cores", n_cores))
} else {
  registerDoSEQ()
  log_message("Sequential processing mode")
}

################################################################################
# SECTION 1: DATA LOADING AND TEMPORAL ALIGNMENT
################################################################################

log_message("Loading aligned multi-platform data...", log_time = TRUE)

# Load datasets
datasets <- list()
data_files <- c(
  "SF" = "SF_UAV_Satellite_aligned_cleaned.csv",
  "EF" = "EF_UAV_Satellite_aligned_cleaned.csv"
)

for (location in names(data_files)) {
  file_path <- file.path(data_dir, data_files[location])
  if (file.exists(file_path)) {
    datasets[[location]] <- read_csv(file_path, show_col_types = FALSE)
    log_message(sprintf("Loaded %s: %d observations", location, nrow(datasets[[location]])))
  } else {
    log_message(sprintf("ERROR: %s not found at %s", location, file_path), "ERROR")
    stop(sprintf("Data file not found: %s", file_path))
  }
}

#' Filter temporal alignment - ONE observation per plot
filter_temporal_alignment <- function(data, max_dap_diff = 5) {
  
  log_message(sprintf("Filtering temporal alignment (max DAP difference: %d days)", max_dap_diff))
  
  # Calculate platform data completeness
  data <- data %>%
    mutate(
      uav_completeness = rowSums(!is.na(select(., starts_with("MS_")))),
      sat_completeness = rowSums(!is.na(select(., starts_with("Satellite_"))))
    )
  
  # Get best observation per plot
  aligned_data <- data %>%
    filter(!is.na(yield)) %>%
    filter(abs(DAP_difference) <= max_dap_diff) %>%
    group_by(Plot_ID) %>%
    arrange(
      abs(DAP_difference),  # First priority: temporal alignment
      desc(uav_completeness + sat_completeness)  # Second: data completeness
    ) %>%
    slice(1) %>%
    ungroup() %>%
    select(-uav_completeness, -sat_completeness)
  
  n_plots_original <- n_distinct(data$Plot_ID)
  n_plots_aligned <- n_distinct(aligned_data$Plot_ID)
  
  log_message(sprintf("  Plots with valid alignment: %d/%d (%.1f%%)", 
                     n_plots_aligned, n_plots_original, 
                     100 * n_plots_aligned / n_plots_original))
  
  return(aligned_data)
}

# Apply temporal filtering to each location
datasets_aligned <- list()
for (name in names(datasets)) {
  log_message(sprintf("\nProcessing %s location:", name))
  data_aligned <- filter_temporal_alignment(datasets[[name]], max_dap_diff = 5)
  datasets_aligned[[name]] <- data_aligned
}

log_message("Data loading and alignment complete", log_time = TRUE)

################################################################################
# SECTION 2: FEATURE EXTRACTION
################################################################################

#' Extract features by platform
extract_features_by_platform <- function(data, platform = "both") {
  
  # Define valid spectral bands and indices patterns
  uav_band_patterns <- c("^MS_Green_", "^MS_Red_", "^MS_Red_edge_", "^MS_NIR_")
  uav_index_patterns <- c("^NDVI_MS_", "^GNDVI_MS_", "^NDRE_MS_", "^NGRDI_MS_",
                          "^SAVI_MS_", "^MSAVI_MS_", "^MTCI_MS_", "^CIgreen_MS_",
                          "^CIrededge_MS_", "^OSAVI_MS_", "^RGRI_MS_")
  
  sat_band_patterns <- c("^Satellite_Satellite_Red_", "^Satellite_Satellite_Green_", 
                         "^Satellite_Satellite_Blue_", "^Satellite_Satellite_NIR_", 
                         "^Satellite_Satellite_Red_edge_", "^Satellite_Satellite_Deep_blue_")
  
  sat_index_patterns <- c("^Satellite_NDVI_Satellite_", "^Satellite_GNDVI_Satellite_", 
                          "^Satellite_NDRE_Satellite_", "^Satellite_GLI_Satellite_", 
                          "^Satellite_NGRDI_Satellite_", "^Satellite_SAVI_Satellite_",
                          "^Satellite_EVI_Satellite_", "^Satellite_MSAVI_Satellite_", 
                          "^Satellite_SIPI_Satellite_", "^Satellite_MTCI_Satellite_", 
                          "^Satellite_CIgreen_Satellite_", "^Satellite_CIrededge_Satellite_",
                          "^Satellite_ARVI_Satellite_", "^Satellite_VARI_Satellite_", 
                          "^Satellite_OSAVI_Satellite_", "^Satellite_TGI_Satellite_", 
                          "^Satellite_ExG_Satellite_", "^Satellite_RGRI_Satellite_")
  
  valid_suffixes <- c("Mean$", "Median$", "P10$", "P25$", "P75$", "P90$")
  
  all_cols <- names(data)
  features <- c()
  
  if (platform == "both" || platform == "combined") {
    patterns <- c(uav_band_patterns, uav_index_patterns, sat_band_patterns, sat_index_patterns)
  } else if (platform == "uav") {
    patterns <- c(uav_band_patterns, uav_index_patterns)
  } else if (platform == "satellite") {
    patterns <- c(sat_band_patterns, sat_index_patterns)
  }
  
  for (pattern in patterns) {
    for (suffix in valid_suffixes) {
      full_pattern <- paste0(pattern, ".*", suffix)
      matches <- all_cols[grepl(full_pattern, all_cols)]
      features <- c(features, matches)
    }
  }
  
  features <- unique(features)
  numeric_features <- features[sapply(data[features], is.numeric)]
  
  metadata_keywords <- c("DAP", "Pass", "Range", "Plot", "Location", "Folder", 
                         "Timepoint", "Growth_Stage", "Acquisition_Date", "Platform",
                         "ProcessedFile", "Spacing", "ID", "LAI_TP", "LAI_DAP", "yield")
  
  for (keyword in metadata_keywords) {
    numeric_features <- numeric_features[!grepl(keyword, numeric_features, ignore.case = TRUE)]
  }
  
  return(numeric_features)
}

extract_platform_features <- function(data) {
  return(list(
    uav = extract_features_by_platform(data, "uav"),
    satellite = extract_features_by_platform(data, "satellite"),
    combined = extract_features_by_platform(data, "both")
  ))
}

################################################################################
# SECTION 3: PLOT-AWARE DATA SPLITTING
################################################################################

#' Create plot-aware train/test split with stratification
create_data_split <- function(data, test_prop = 0.3, seed = 123) {
  
  set.seed(seed)
  
  # Verify one observation per plot
  obs_per_plot <- data %>% count(Plot_ID) %>% pull(n)
  if (any(obs_per_plot > 1)) {
    stop("Data contains multiple observations per plot! Apply temporal filtering first.")
  }
  
  # Get unique plots with metadata for stratification
  plot_info <- data %>%
    select(Plot_ID, Spacing, Leaf_type, yield) %>%
    distinct()
  
  # Create stratification bins
  n_unique_yield <- length(unique(plot_info$yield))
  n_bins <- min(5, n_unique_yield)
  
  if (n_bins > 1) {
    plot_info <- plot_info %>%
      mutate(
        yield_bin = cut(yield, breaks = n_bins, labels = FALSE, include.lowest = TRUE),
        strata = paste(Spacing, Leaf_type, yield_bin, sep = "_")
      )
  } else {
    plot_info <- plot_info %>%
      mutate(
        yield_bin = 1,
        strata = paste(Spacing, Leaf_type, sep = "_")
      )
  }
  
  # Stratified split
  test_plots <- plot_info %>%
    group_by(strata) %>%
    mutate(n_strata = n()) %>%
    ungroup() %>%
    group_by(strata) %>%
    sample_n(size = pmax(1, round(n() * test_prop))) %>%
    ungroup() %>%
    pull(Plot_ID)
  
  train_plots <- setdiff(plot_info$Plot_ID, test_plots)
  
  # Create datasets
  train_data <- data %>% filter(Plot_ID %in% train_plots)
  test_data <- data %>% filter(Plot_ID %in% test_plots)
  
  log_message(sprintf("  Plot-aware split - Train: %d plots, Test: %d plots",
                     length(train_plots), length(test_plots)))
  
  # Verify no plot overlap
  if (length(intersect(train_plots, test_plots)) > 0) {
    stop("Plot overlap detected in train/test split!")
  }
  
  return(list(
    train = train_data, 
    test = test_data,
    train_plots = train_plots,
    test_plots = test_plots
  ))
}

################################################################################
# SECTION 4: SAMPLING FUNCTIONS
################################################################################

#' Plot-aware stratified sampling for bootstrap
stratified_sample_plots <- function(data, prop = 1.0, target_var = "yield", n_strata = 5) {
  
  plot_data <- data %>%
    group_by(Plot_ID) %>%
    summarise(
      !!target_var := first(!!sym(target_var)),
      Spacing = first(Spacing),
      Leaf_type = first(Leaf_type),
      .groups = "drop"
    )
  
  plot_data <- plot_data %>%
    mutate(
      target_quantile = cut(!!sym(target_var), 
                           breaks = quantile(!!sym(target_var), 
                                           probs = seq(0, 1, 1/n_strata),
                                           na.rm = TRUE),
                           labels = FALSE, include.lowest = TRUE)
    )
  
  sampled_plots <- plot_data %>%
    group_by(target_quantile, Spacing) %>%
    sample_frac(size = prop, replace = TRUE) %>%
    pull(Plot_ID)
  
  plot_counts <- table(sampled_plots)
  
  sampled_indices <- c()
  for (plot in names(plot_counts)) {
    plot_idx <- which(data$Plot_ID == plot)
    sampled_indices <- c(sampled_indices, rep(plot_idx, plot_counts[plot]))
  }
  
  return(sampled_indices)
}

################################################################################
# SECTION 5: EVALUATION FUNCTIONS
################################################################################

#' Calculate enhanced metrics
calculate_metrics_enhanced <- function(observed, predicted) {
  
  valid_idx <- !is.na(observed) & !is.na(predicted)
  obs <- observed[valid_idx]
  pred <- predicted[valid_idx]
  
  if (length(obs) < 5) {
    return(list(R2 = NA, RMSE = NA, MAE = NA, Bias = NA, MAPE = NA, n = length(obs)))
  }
  
  residuals <- obs - pred
  ss_res <- sum(residuals^2)
  ss_tot <- sum((obs - mean(obs))^2)
  
  metrics <- list(
    R2 = 1 - ss_res / ss_tot,
    RMSE = sqrt(mean(residuals^2)),
    MAE = mean(abs(residuals)),
    Bias = mean(residuals),
    MAPE = mean(abs(residuals / (obs + 0.01))) * 100,
    n = length(obs)
  )
  
  return(metrics)
}

################################################################################
# SECTION 6: COMPLEMENTARITY ANALYSIS
################################################################################

#' Enhanced complementarity analysis with TEST SET evaluation
analyze_complementarity <- function(train_data, test_data, features, n_bootstrap = N_BOOTSTRAP) {
  
  log_message("\n=== Platform Complementarity Analysis ===")
  log_message(sprintf("  Running %d bootstrap iterations...", n_bootstrap), log_time = TRUE)
  log_message("  All models evaluated on test set for fair comparison")
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
  
  # Bootstrap results storage
  bootstrap_results <- foreach(b = 1:n_bootstrap,
                              .packages = c("tidyverse", "caret", "randomForest"),
                              .export = c("stratified_sample_plots", "extract_features_by_platform",
                                         "extract_platform_features", "features", 
                                         "nearZeroVar", "preProcess"),
                              .combine = 'rbind',
                              .errorhandling = 'remove') %dopar% {
    
    set.seed(b * 100)
    
    # Force garbage collection
    if (b %% 10 == 0) {
      gc(verbose = FALSE)
    }
    
    # Plot-aware bootstrap sample of training data
    boot_idx <- stratified_sample_plots(train_data, 1.0, "yield")
    train_boot <- train_data[boot_idx, ]
    
    # Prepare TRAINING data
    X_uav <- train_boot %>% select(all_of(features$uav)) %>% select_if(is.numeric)
    X_sat <- train_boot %>% select(all_of(features$satellite)) %>% select_if(is.numeric)
    X_combined <- train_boot %>% select(all_of(features$combined)) %>% select_if(is.numeric)
    y <- train_boot$yield
    
    # Prepare TEST data
    X_test_uav <- test_data %>% select(all_of(names(X_uav))) %>% select_if(is.numeric)
    X_test_sat <- test_data %>% select(all_of(names(X_sat))) %>% select_if(is.numeric)
    X_test_combined <- test_data %>% select(all_of(names(X_combined))) %>% select_if(is.numeric)
    y_test <- test_data$yield
    
    if (ncol(X_uav) < 5 || ncol(X_sat) < 5) return(NULL)
    
    # Remove near-zero variance features
    nzv_uav <- nearZeroVar(X_uav)
    nzv_sat <- nearZeroVar(X_sat)
    nzv_combined <- nearZeroVar(X_combined)
    
    if (length(nzv_uav) > 0) {
      X_uav <- X_uav[, -nzv_uav]
      X_test_uav <- X_test_uav[, -nzv_uav]
    }
    if (length(nzv_sat) > 0) {
      X_sat <- X_sat[, -nzv_sat]
      X_test_sat <- X_test_sat[, -nzv_sat]
    }
    if (length(nzv_combined) > 0) {
      X_combined <- X_combined[, -nzv_combined]
      X_test_combined <- X_test_combined[, -nzv_combined]
    }
    
    # Preprocess training data
    preproc_uav <- preProcess(X_uav, method = c("center", "scale", "medianImpute"))
    X_uav_proc <- predict(preproc_uav, X_uav)
    
    preproc_sat <- preProcess(X_sat, method = c("center", "scale", "medianImpute"))
    X_sat_proc <- predict(preproc_sat, X_sat)
    
    preproc_combined <- preProcess(X_combined, method = c("center", "scale", "medianImpute"))
    X_combined_proc <- predict(preproc_combined, X_combined)
    
    # Preprocess test data with same transformations
    X_test_uav_proc <- predict(preproc_uav, X_test_uav)
    X_test_sat_proc <- predict(preproc_sat, X_test_sat)
    X_test_combined_proc <- predict(preproc_combined, X_test_combined)
    
    # Train models
    rf_uav <- randomForest(X_uav_proc, y, ntree = 300, importance = TRUE)
    rf_sat <- randomForest(X_sat_proc, y, ntree = 300, importance = TRUE)
    rf_combined <- randomForest(X_combined_proc, y, ntree = 300, importance = TRUE)
    
    # Calculate OOB metrics (for reference only)
    r2_uav_oob <- 1 - mean((rf_uav$predicted - y)^2) / var(y)
    r2_sat_oob <- 1 - mean((rf_sat$predicted - y)^2) / var(y)
    r2_combined_oob <- 1 - mean((rf_combined$predicted - y)^2) / var(y)
    
    # Calculate TEST SET metrics for fair comparison
    pred_uav_test <- predict(rf_uav, X_test_uav_proc)
    pred_sat_test <- predict(rf_sat, X_test_sat_proc)
    pred_combined_test <- predict(rf_combined, X_test_combined_proc)
    
    # TEST SET R²
    r2_uav_test <- 1 - mean((pred_uav_test - y_test)^2) / var(y_test)
    r2_sat_test <- 1 - mean((pred_sat_test - y_test)^2) / var(y_test)
    r2_combined_test <- 1 - mean((pred_combined_test - y_test)^2) / var(y_test)
    
    # Complementarity on TEST SET
    complementarity_score <- r2_combined_test - max(r2_uav_test, r2_sat_test)
    synergy_ratio <- r2_combined_test / (r2_uav_test + r2_sat_test + 1e-10)
    
    # Residual correlation on TEST SET
    residuals_uav_test <- y_test - pred_uav_test
    residuals_sat_test <- y_test - pred_sat_test
    residual_correlation <- cor(residuals_uav_test, residuals_sat_test)
    
    # Store feature importance (top 20 features)
    uav_importance <- importance(rf_uav)[, "%IncMSE"]
    sat_importance <- importance(rf_sat)[, "%IncMSE"]
    combined_importance <- importance(rf_combined)[, "%IncMSE"]
    
    data.frame(
      bootstrap = b,
      # OOB metrics (for reference)
      r2_uav_oob = r2_uav_oob,
      r2_sat_oob = r2_sat_oob,
      r2_combined_oob = r2_combined_oob,
      # TEST SET metrics (for fair comparison)
      r2_uav = r2_uav_test,
      r2_sat = r2_sat_test,
      r2_combined = r2_combined_test,
      complementarity_score = complementarity_score,
      synergy_ratio = synergy_ratio,
      residual_correlation = residual_correlation,
      n_plots = n_distinct(train_boot$Plot_ID),
      n_features_uav = ncol(X_uav_proc),
      n_features_sat = ncol(X_sat_proc),
      n_features_combined = ncol(X_combined_proc)
    )
  }
  
  close(pb)
  setTxtProgressBar(pb, n_bootstrap)
  
  log_message("  Complementarity analysis complete", log_time = TRUE)
  
  # Calculate statistics with CI
  results <- list()
  
  # TEST SET metrics (primary results for fair comparison)
  for (metric in c("r2_uav", "r2_sat", "r2_combined", "complementarity_score", 
                   "synergy_ratio", "residual_correlation")) {
    values <- bootstrap_results[[metric]]
    results[[metric]] <- list(
      mean = mean(values, na.rm = TRUE),
      sd = sd(values, na.rm = TRUE),
      ci_lower = quantile(values, (1 - CONFIDENCE_LEVEL) / 2, na.rm = TRUE),
      ci_upper = quantile(values, 1 - (1 - CONFIDENCE_LEVEL) / 2, na.rm = TRUE),
      n = sum(!is.na(values))
    )
  }
  
  # OOB metrics (for comparison to show bias)
  for (metric in c("r2_uav_oob", "r2_sat_oob", "r2_combined_oob")) {
    values <- bootstrap_results[[metric]]
    results[[metric]] <- list(
      mean = mean(values, na.rm = TRUE),
      sd = sd(values, na.rm = TRUE),
      ci_lower = quantile(values, (1 - CONFIDENCE_LEVEL) / 2, na.rm = TRUE),
      ci_upper = quantile(values, 1 - (1 - CONFIDENCE_LEVEL) / 2, na.rm = TRUE),
      n = sum(!is.na(values))
    )
  }
  
  results$bootstrap_results <- bootstrap_results
  results$mean_plots_per_bootstrap <- mean(bootstrap_results$n_plots)
  
  log_message(sprintf("  TEST SET Complementarity score: %.3f [%.3f, %.3f]",
                     results$complementarity_score$mean,
                     results$complementarity_score$ci_lower,
                     results$complementarity_score$ci_upper))
  
  # Log the difference between OOB and test set
  log_message(sprintf("  OOB vs Test R² difference: UAV=%.3f, Sat=%.3f, Combined=%.3f",
                     results$r2_uav_oob$mean - results$r2_uav$mean,
                     results$r2_sat_oob$mean - results$r2_sat$mean,
                     results$r2_combined_oob$mean - results$r2_combined$mean))
  
  return(results)
}

################################################################################
# SECTION 7: FEATURE IMPORTANCE ANALYSIS
################################################################################

#' Analyze feature importance across platforms
analyze_feature_importance <- function(train_data, test_data, features, n_bootstrap = 50) {
  
  log_message("\n=== Feature Importance Analysis ===")
  log_message(sprintf("  Running %d bootstrap iterations...", n_bootstrap))
  
  # Store importance for each bootstrap
  importance_results <- list(
    uav = list(),
    satellite = list(),
    combined = list()
  )
  
  for (b in 1:n_bootstrap) {
    set.seed(b * 200)
    
    # Bootstrap sample
    boot_idx <- stratified_sample_plots(train_data, 1.0, "yield")
    train_boot <- train_data[boot_idx, ]
    
    # Train models for each platform
    for (platform in c("uav", "satellite", "combined")) {
      X <- train_boot %>% select(all_of(features[[platform]])) %>% select_if(is.numeric)
      y <- train_boot$yield
      
      if (ncol(X) < 5) next
      
      # Remove near-zero variance
      nzv <- nearZeroVar(X)
      if (length(nzv) > 0) X <- X[, -nzv]
      
      # Preprocess
      preproc <- preProcess(X, method = c("center", "scale", "medianImpute"))
      X_proc <- predict(preproc, X)
      
      # Train model
      rf_model <- randomForest(X_proc, y, ntree = 300, importance = TRUE)
      
      # Store importance
      imp_matrix <- importance(rf_model)
      importance_results[[platform]][[b]] <- data.frame(
        feature = rownames(imp_matrix),
        IncMSE = imp_matrix[, "%IncMSE"],
        IncNodePurity = imp_matrix[, "IncNodePurity"],
        bootstrap = b
      )
    }
  }
  
  # Aggregate importance
  aggregated_importance <- list()
  
  for (platform in names(importance_results)) {
    if (length(importance_results[[platform]]) > 0) {
      # Combine all bootstraps
      all_imp <- do.call(rbind, importance_results[[platform]])
      
      # Calculate statistics
      aggregated_importance[[platform]] <- all_imp %>%
        group_by(feature) %>%
        summarise(
          IncMSE_mean = mean(IncMSE),
          IncMSE_sd = sd(IncMSE),
          IncMSE_ci_lower = quantile(IncMSE, (1 - CONFIDENCE_LEVEL) / 2),
          IncMSE_ci_upper = quantile(IncMSE, 1 - (1 - CONFIDENCE_LEVEL) / 2),
          IncNodePurity_mean = mean(IncNodePurity),
          IncNodePurity_sd = sd(IncNodePurity),
          n_bootstraps = n(),
          .groups = "drop"
        ) %>%
        arrange(desc(IncMSE_mean))
    }
  }
  
  log_message("  Feature importance analysis complete")
  
  return(aggregated_importance)
}

################################################################################
# SECTION 8: STATISTICAL TESTING
################################################################################

#' Perform statistical comparisons between platforms
compare_platforms_statistically <- function(comp_results) {
  
  bootstrap_data <- comp_results$bootstrap_results
  
  # Wilcoxon signed-rank tests
  tests <- list()
  
  # UAV vs Satellite
  test_uav_sat <- wilcox.test(bootstrap_data$r2_uav, bootstrap_data$r2_sat, 
                             paired = TRUE, exact = FALSE)
  tests$uav_vs_sat <- list(
    p_value = test_uav_sat$p.value,
    statistic = test_uav_sat$statistic,
    effect_size = median(bootstrap_data$r2_uav - bootstrap_data$r2_sat)
  )
  
  # Combined vs best single platform
  best_single <- pmax(bootstrap_data$r2_uav, bootstrap_data$r2_sat)
  test_combined <- wilcox.test(bootstrap_data$r2_combined, best_single, 
                              paired = TRUE, exact = FALSE)
  tests$combined_vs_best <- list(
    p_value = test_combined$p.value,
    statistic = test_combined$statistic,
    effect_size = median(bootstrap_data$r2_combined - best_single)
  )
  
  # Test if complementarity is significantly positive
  test_comp <- wilcox.test(bootstrap_data$complementarity_score, mu = 0, 
                          alternative = "greater", exact = FALSE)
  tests$complementarity_positive <- list(
    p_value = test_comp$p.value,
    statistic = test_comp$statistic,
    effect_size = median(bootstrap_data$complementarity_score)
  )
  
  return(tests)
}

################################################################################
# SECTION 9: VISUALIZATION FUNCTIONS
################################################################################

#' Plot complementarity results
plot_complementarity_results <- function(comp_results, stat_tests, location, output_dir) {
  
  pdf(file.path(output_dir, paste0(location, "_complementarity_yield_analysis.pdf")), 
      width = 16, height = 12)
  
  # 1. Platform comparison with TEST SET performance
  platform_data <- data.frame(
    Platform = c("UAV Only", "Satellite Only", "Combined"),
    R2_mean = c(comp_results$r2_uav$mean, 
                comp_results$r2_sat$mean, 
                comp_results$r2_combined$mean),
    R2_lower = c(comp_results$r2_uav$ci_lower,
                 comp_results$r2_sat$ci_lower,
                 comp_results$r2_combined$ci_lower),
    R2_upper = c(comp_results$r2_uav$ci_upper,
                 comp_results$r2_sat$ci_upper,
                 comp_results$r2_combined$ci_upper)
  )
  
  # Add significance annotations
  sig_label <- function(p) {
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("ns")
  }
  
  p1 <- ggplot(platform_data, aes(x = Platform, y = R2_mean, fill = Platform)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = R2_lower, ymax = R2_upper), width = 0.2, size = 1) +
    geom_text(aes(label = sprintf("%.3f", R2_mean)), vjust = -1.5) +
    theme_minimal(base_size = 12) +
    labs(title = paste0(location, " - Test Set Performance by Platform"),
         subtitle = sprintf("Based on %d bootstrap iterations", comp_results$r2_uav$n),
         y = "Test Set R²") +
    scale_fill_manual(values = c("UAV Only" = "#1f77b4", 
                                 "Satellite Only" = "#ff7f0e",
                                 "Combined" = "#2ca02c")) +
    ylim(0, max(platform_data$R2_upper) * 1.2) +
    theme(legend.position = "none") +
    annotate("text", x = 1.5, y = max(platform_data$R2_upper) * 1.15, 
             label = sprintf("p = %.3f %s", stat_tests$uav_vs_sat$p_value,
                            sig_label(stat_tests$uav_vs_sat$p_value)))
  
  # 2. Bootstrap distributions
  boot_data <- comp_results$bootstrap_results
  boot_long <- boot_data %>%
    select(bootstrap, r2_uav, r2_sat, r2_combined) %>%
    pivot_longer(cols = -bootstrap, names_to = "Platform", values_to = "R2") %>%
    mutate(Platform = recode(Platform,
                             r2_uav = "UAV Only",
                             r2_sat = "Satellite Only",
                             r2_combined = "Combined"))
  
  p2 <- ggplot(boot_long, aes(x = R2, fill = Platform)) +
    geom_density(alpha = 0.6) +
    facet_wrap(~Platform, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(title = "Bootstrap Distribution of Test Set R² Values",
         x = "Test Set R²", y = "Density") +
    scale_fill_manual(values = c("UAV Only" = "#1f77b4", 
                                 "Satellite Only" = "#ff7f0e",
                                 "Combined" = "#2ca02c")) +
    theme(legend.position = "none")
  
  # 3. Complementarity metrics
  comp_metrics <- data.frame(
    Metric = c("Complementarity\nScore", "Synergy\nRatio", "Residual\nCorrelation"),
    Value_mean = c(comp_results$complementarity_score$mean, 
                   comp_results$synergy_ratio$mean,
                   comp_results$residual_correlation$mean),
    Value_lower = c(comp_results$complementarity_score$ci_lower,
                    comp_results$synergy_ratio$ci_lower,
                    comp_results$residual_correlation$ci_lower),
    Value_upper = c(comp_results$complementarity_score$ci_upper,
                    comp_results$synergy_ratio$ci_upper,
                    comp_results$residual_correlation$ci_upper)
  )
  
  p3 <- ggplot(comp_metrics, aes(x = Metric, y = Value_mean, fill = Metric)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = Value_lower, ymax = Value_upper), 
                  width = 0.2, linewidth = 1) +
    geom_text(aes(label = sprintf("%.3f", Value_mean)), vjust = -1) +
    theme_minimal(base_size = 12) +
    labs(title = "Platform Complementarity Metrics",
         subtitle = sprintf("Complementarity test: p = %.3f %s", 
                           stat_tests$complementarity_positive$p_value,
                           sig_label(stat_tests$complementarity_positive$p_value)),
         y = "Value") +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Set2")
  
  # 4. OOB vs Test Set comparison
  comparison_data <- data.frame(
    Platform = rep(c("UAV", "Satellite", "Combined"), 2),
    Evaluation = c(rep("OOB", 3), rep("Test Set", 3)),
    R2 = c(comp_results$r2_uav_oob$mean, comp_results$r2_sat_oob$mean, 
           comp_results$r2_combined_oob$mean,
           comp_results$r2_uav$mean, comp_results$r2_sat$mean, 
           comp_results$r2_combined$mean)
  )
  
  p4 <- ggplot(comparison_data, aes(x = Platform, y = R2, fill = Evaluation)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    theme_minimal() +
    labs(title = "OOB vs Test Set Performance",
         subtitle = "Showing optimistic bias in OOB estimates",
         y = "R²") +
    scale_fill_manual(values = c("OOB" = "#fdae61", "Test Set" = "#2c7bb6"))
  
  # Combine plots
  print((p1 | p2) / (p3 | p4))
  
  dev.off()
}

#' Plot feature importance
plot_feature_importance <- function(importance_results, location, output_dir) {
  
  pdf(file.path(output_dir, paste0(location, "_feature_importance_yield.pdf")), 
      width = 14, height = 10)
  
  for (platform in names(importance_results)) {
    if (!is.null(importance_results[[platform]]) && 
        nrow(importance_results[[platform]]) > 0) {
      
      # Get top 20 features
      top_features <- importance_results[[platform]] %>%
        slice_head(n = 20)
      
      p <- ggplot(top_features, 
                  aes(x = reorder(feature, IncMSE_mean), y = IncMSE_mean)) +
        geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.8) +
        geom_errorbar(aes(ymin = IncMSE_ci_lower, ymax = IncMSE_ci_upper),
                      width = 0.2) +
        coord_flip() +
        theme_minimal(base_size = 10) +
        labs(title = paste0("Top 20 Features - ", toupper(platform), " Platform"),
             subtitle = paste0(location, " - With 95% confidence intervals"),
             x = "", y = "Mean %IncMSE") +
        theme(axis.text.y = element_text(size = 8))
      
      print(p)
    }
  }
  
  dev.off()
}

################################################################################
# SECTION 10: SAVE RAW OUTPUTS
################################################################################

#' Save all raw outputs
save_raw_outputs <- function(results, location, output_dir) {
  
  log_message(sprintf("  Saving raw output files for %s...", location))
  
  # 1. Save bootstrap results
  if (!is.null(results$complementarity)) {
    comp_boot_file <- file.path(output_dir, "raw_outputs", location,
                                "bootstrap_results", "complementarity_bootstrap_results.csv")
    write.csv(results$complementarity$bootstrap_results, comp_boot_file, row.names = FALSE)
    
    # Platform comparison summary
    platform_summary <- data.frame(
      Platform = c("UAV_Only", "Satellite_Only", "Combined"),
      R2_mean = c(results$complementarity$r2_uav$mean,
                  results$complementarity$r2_sat$mean,
                  results$complementarity$r2_combined$mean),
      R2_sd = c(results$complementarity$r2_uav$sd,
                results$complementarity$r2_sat$sd,
                results$complementarity$r2_combined$sd),
      R2_ci_lower = c(results$complementarity$r2_uav$ci_lower,
                      results$complementarity$r2_sat$ci_lower,
                      results$complementarity$r2_combined$ci_lower),
      R2_ci_upper = c(results$complementarity$r2_uav$ci_upper,
                      results$complementarity$r2_sat$ci_upper,
                      results$complementarity$r2_combined$ci_upper)
    )
    
    platform_file <- file.path(output_dir, "raw_outputs", location,
                               "platform_comparison", "platform_comparison_summary.csv")
    write.csv(platform_summary, platform_file, row.names = FALSE)
  }
  
  # 2. Save feature importance
  if (!is.null(results$feature_importance)) {
    for (platform in names(results$feature_importance)) {
      if (!is.null(results$feature_importance[[platform]])) {
        imp_file <- file.path(output_dir, "raw_outputs", location,
                              "feature_importance", 
                              paste0(platform, "_feature_importance.csv"))
        write.csv(results$feature_importance[[platform]], imp_file, row.names = FALSE)
      }
    }
  }
  
  # 3. Save statistical test results
  if (!is.null(results$statistical_tests)) {
    stat_tests <- results$statistical_tests
    stat_results <- data.frame(
      Test = c("UAV_vs_Satellite", "Combined_vs_Best_Single", "Complementarity_Positive"),
      p_value = c(stat_tests$uav_vs_sat$p_value,
                  stat_tests$combined_vs_best$p_value,
                  stat_tests$complementarity_positive$p_value),
      test_statistic = c(stat_tests$uav_vs_sat$statistic,
                         stat_tests$combined_vs_best$statistic,
                         stat_tests$complementarity_positive$statistic),
      effect_size = c(stat_tests$uav_vs_sat$effect_size,
                      stat_tests$combined_vs_best$effect_size,
                      stat_tests$complementarity_positive$effect_size)
    )
    
    stat_file <- file.path(output_dir, "summary_tables",
                           paste0(location, "_statistical_test_results.csv"))
    write.csv(stat_results, stat_file, row.names = FALSE)
  }
  
  log_message(sprintf("  Raw output files saved for %s", location))
}

################################################################################
# SECTION 11: MAIN EXECUTION
################################################################################

log_message("\n=== MAIN EXECUTION ===", log_time = TRUE)

# Process each location
all_results <- list()
all_timing <- list()

for (location in names(datasets_aligned)) {
  
  location_start_time <- Sys.time()
  log_message(sprintf("\n=== Processing %s Location ===", location), log_time = TRUE)
  
  data <- datasets_aligned[[location]]
  
  # Extract features
  features <- extract_platform_features(data)
  log_message(sprintf("  Features - UAV: %d, Satellite: %d, Combined: %d",
                     length(features$uav), length(features$satellite), 
                     length(features$combined)))
  
  # Create plot-aware train/test split
  splits <- create_data_split(data, test_prop = 0.3)
  
  # Initialize results
  results <- list()
  
  # Run complementarity analysis
  results$complementarity <- analyze_complementarity(
    splits$train, splits$test, features, n_bootstrap = N_BOOTSTRAP
  )
  
  # Run feature importance analysis
  results$feature_importance <- analyze_feature_importance(
    splits$train, splits$test, features, n_bootstrap = 50
  )
  
  # Statistical tests
  results$statistical_tests <- compare_platforms_statistically(results$complementarity)
  
  # Store results
  all_results[[location]] <- results
  
  # Create visualizations
  viz_dir <- file.path(phase5_dir, "visualizations")
  plot_complementarity_results(results$complementarity, results$statistical_tests, 
                               location, viz_dir)
  plot_feature_importance(results$feature_importance, location, viz_dir)
  
  # Save raw outputs
  save_raw_outputs(results, location, phase5_dir)
  
  # Save complete results
  saveRDS(results, file.path(phase5_dir, "model_performance", 
                            sprintf("%s_complementarity_results_yield.rds", location)))
  
  # Calculate and log location processing time
  location_time <- as.numeric(difftime(Sys.time(), location_start_time, units = "mins"))
  all_timing[[location]] <- location_time
  log_message(sprintf("  %s processing complete. Time: %.1f minutes", location, location_time), 
              log_time = TRUE)
}

################################################################################
# SECTION 12: CREATE SUMMARY REPORT
################################################################################

log_message("\n=== Creating Summary Report ===", log_time = TRUE)

# Create comprehensive summary
summary_data <- data.frame()

for (location in names(all_results)) {
  results <- all_results[[location]]
  
  summary_data <- rbind(summary_data, data.frame(
    Location = location,
    UAV_R2 = results$complementarity$r2_uav$mean,
    UAV_R2_CI = sprintf("[%.3f, %.3f]", 
                        results$complementarity$r2_uav$ci_lower,
                        results$complementarity$r2_uav$ci_upper),
    Satellite_R2 = results$complementarity$r2_sat$mean,
    Satellite_R2_CI = sprintf("[%.3f, %.3f]", 
                              results$complementarity$r2_sat$ci_lower,
                              results$complementarity$r2_sat$ci_upper),
    Combined_R2 = results$complementarity$r2_combined$mean,
    Combined_R2_CI = sprintf("[%.3f, %.3f]", 
                             results$complementarity$r2_combined$ci_lower,
                             results$complementarity$r2_combined$ci_upper),
    Complementarity_Score = results$complementarity$complementarity_score$mean,
    Complementarity_CI = sprintf("[%.3f, %.3f]", 
                                 results$complementarity$complementarity_score$ci_lower,
                                 results$complementarity$complementarity_score$ci_upper),
    stringsAsFactors = FALSE
  ))
}

# Save summary
write_csv(summary_data, file.path(phase5_dir, "summary_tables", 
                                  "platform_comparison_summary_yield.csv"))

# Save timing information
timing_df <- data.frame(
  Location = names(all_timing),
  Processing_time_mins = unlist(all_timing),
  Total_time_mins = as.numeric(difftime(Sys.time(), overall_start_time, units = "mins"))
)
write_csv(timing_df, file.path(phase5_dir, "time_logs", "processing_times_summary.csv"))

# Create text report
sink(file.path(phase5_dir, "complementarity_analysis_summary_report_yield.txt"))

cat("MULTI-PLATFORM COMPLEMENTARITY ANALYSIS FOR YIELD - SUMMARY\n")
cat("===========================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total execution time:", sprintf("%.1f minutes\n", 
    as.numeric(difftime(Sys.time(), overall_start_time, units = "mins"))))
cat("\n")

cat("ANALYSIS CONFIGURATION:\n")
cat(sprintf("- Bootstrap iterations: %d\n", N_BOOTSTRAP))
cat(sprintf("- Confidence level: %.0f%%\n", CONFIDENCE_LEVEL * 100))
cat("- Temporal alignment: Max 5 days UAV-Satellite difference\n")
cat("- Plot-aware: One observation per plot, no data leakage\n")
cat("- Evaluation: All models tested on same held-out test set\n")
cat(sprintf("- Parallel processing: %s with %d cores\n", 
            ifelse(USE_PARALLEL, "Enabled", "Disabled"),
            ifelse(USE_PARALLEL, n_cores, 1)))
cat("\n")

# Results by location
for (location in names(all_results)) {
  results <- all_results[[location]]
  
  cat(sprintf("\n%s LOCATION:\n", toupper(location)))
  cat(sprintf("  Processing time: %.1f minutes\n", all_timing[[location]]))
  
  cat("\n  Platform Performance (Test Set):\n")
  cat(sprintf("    UAV only R²: %.3f [%.3f, %.3f]\n",
              results$complementarity$r2_uav$mean,
              results$complementarity$r2_uav$ci_lower,
              results$complementarity$r2_uav$ci_upper))
  cat(sprintf("    Satellite only R²: %.3f [%.3f, %.3f]\n",
              results$complementarity$r2_sat$mean,
              results$complementarity$r2_sat$ci_lower,
              results$complementarity$r2_sat$ci_upper))
  cat(sprintf("    Combined R²: %.3f [%.3f, %.3f]\n",
              results$complementarity$r2_combined$mean,
              results$complementarity$r2_combined$ci_lower,
              results$complementarity$r2_combined$ci_upper))
  
  cat("\n  Complementarity Metrics:\n")
  cat(sprintf("    Complementarity score: %.3f [%.3f, %.3f]\n",
              results$complementarity$complementarity_score$mean,
              results$complementarity$complementarity_score$ci_lower,
              results$complementarity$complementarity_score$ci_upper))
  cat(sprintf("    Synergy ratio: %.3f\n", 
              results$complementarity$synergy_ratio$mean))
  cat(sprintf("    Residual correlation: %.3f\n", 
              results$complementarity$residual_correlation$mean))
  
  cat("\n  Statistical Tests:\n")
  cat(sprintf("    UAV vs Satellite: p = %.4f\n", 
              results$statistical_tests$uav_vs_sat$p_value))
  cat(sprintf("    Combined vs Best Single: p = %.4f\n", 
              results$statistical_tests$combined_vs_best$p_value))
  cat(sprintf("    Complementarity > 0: p = %.4f\n", 
              results$statistical_tests$complementarity_positive$p_value))
}

cat("\n\nFILES GENERATED:\n")
cat(sprintf("- Main results: %s\n", phase5_dir))
cat("- Visualizations: visualizations/\n")
cat("- Raw outputs: raw_outputs/[LOCATION]/\n")
cat("- Summary tables: summary_tables/\n")

sink()

log_message("\n=== SIMPLIFIED MULTI-PLATFORM COMPLEMENTARITY ANALYSIS FOR YIELD COMPLETE ===", log_time = TRUE)
log_message(sprintf("Total execution time: %.1f minutes", 
                   as.numeric(difftime(Sys.time(), overall_start_time, units = "mins"))))
log_message(sprintf("All results saved to: %s", phase5_dir))

# Cleanup
if (USE_PARALLEL && exists("cl")) {
  stopCluster(cl)
  log_message("Parallel cluster closed")
}

# End of script