################################################################################
# Phase 4: Growth Stage Analysis for Peak LAI Prediction - TEMPORAL DATA VERSION
# 
# Purpose: Analyze peak LAI prediction performance across growth stages
#          using raw temporal data to avoid data leakage
#
# Key Changes:
# - Uses raw temporal data files (*_with_DAP_cleaned.csv)
# - Dynamically aggregates features within growth stages
# - Removes LAI trajectory features that contain target information
# - Maintains all statistical rigor (bootstrap, repeated holdout, CI)
#
# Growth Stage Definitions:
# - Based on Growth_Stage column in temporal data
# - Features aggregated within each stage dynamically
#
# Models: lasso, ridge, lm, pls, rf, xgboost, svm
#
# Author: Soybean Remote Sensing Analysis Pipeline - Temporal Version
# Date: 2024
################################################################################

# Load required libraries
cat("Loading required libraries...\n")

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(glmnet)
  library(randomForest)
  library(xgboost)
  library(e1071)
  library(pls)
  library(corrplot)
  library(ggplot2)
  library(viridis)
  library(gridExtra)
  library(ggrepel)
  library(reshape2)
  library(patchwork)
  library(parallel)
  library(doParallel)
  library(boot)
})

# Set seed and options
set.seed(123)
options(scipen = 999)

# Configuration
DEBUG_MODE <- FALSE  # Set to TRUE for quick debugging
USE_PARALLEL <- TRUE  # Set to FALSE if parallel issues on cluster
N_BOOTSTRAP <- ifelse(DEBUG_MODE, 10, 100)  # Number of bootstrap iterations
N_REPETITIONS <- ifelse(DEBUG_MODE, 5, 50)  # Number of repeated holdout splits
CONFIDENCE_LEVEL <- 0.95  # For confidence intervals

# Check if on cluster
on_cluster <- Sys.getenv("SLURM_JOB_ID") != ""
if (on_cluster) {
  cat(sprintf("Running on cluster - Job ID: %s\n", Sys.getenv("SLURM_JOB_ID")))
}

################################################################################
# SECTION 0: SETUP AND PATHS
################################################################################

# Set paths
current_dir <- getwd()
project_root <- ifelse(basename(current_dir) == "scripts", dirname(current_dir), current_dir)

# Direct paths
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "soybean_analysis_outputs", "phase4_growth_stage_peak_lai_temporal")

# Create output directories
subdirs <- c("stage_specific_models", "feature_importance", 
             "temporal_accumulation", "platform_comparison", 
             "stage_combinations", "bootstrap_results", 
             "confidence_intervals", "statistical_tests", "tables")
for (subdir in subdirs) {
  dir_path <- file.path(output_dir, subdir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

# Logging function
log_file <- file.path(output_dir, "phase4_growth_stage_peak_lai_temporal_log.txt")
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- sprintf("[%s] %s: %s", timestamp, level, msg)
  cat(formatted_msg, "\n")
  cat(formatted_msg, "\n", file = log_file, append = TRUE)
}

log_message("Starting Phase 4: Growth Stage Analysis for Peak LAI - TEMPORAL VERSION")
log_message(sprintf("Configuration: %d bootstrap iterations, %d repeated holdouts", 
                    N_BOOTSTRAP, N_REPETITIONS))

################################################################################
# SECTION 1: LOAD TEMPORAL DATA
################################################################################

log_message("Loading temporal data with DAP information...")

# Define temporal data files
temporal_files <- c(
  "SF_UAV" = "SF_UAV_with_DAP_cleaned.csv",
  "EF_UAV" = "EF_UAV_with_DAP_cleaned.csv",
  "SF_Satellite" = "SF_Satellite_with_DAP_cleaned.csv",
  "EF_Satellite" = "EF_Satellite_with_DAP_cleaned.csv"
)

# Load temporal datasets
temporal_datasets <- list()
for (dataset_name in names(temporal_files)) {
  file_path <- file.path(data_dir, temporal_files[dataset_name])
  if (file.exists(file_path)) {
    temporal_datasets[[dataset_name]] <- read_csv(file_path, show_col_types = FALSE)
    log_message(sprintf("Loaded %s: %d observations", dataset_name, 
                        nrow(temporal_datasets[[dataset_name]])))
    
    # Check growth stages
    stages <- unique(temporal_datasets[[dataset_name]]$Growth_Stage)
    log_message(sprintf("  Growth stages in %s: %s", dataset_name, 
                        paste(stages, collapse = ", ")))
  } else {
    log_message(sprintf("Warning: %s not found", file_path), "WARNING")
  }
}

################################################################################
# SECTION 2: TEMPORAL DATA PREPARATION FUNCTIONS
################################################################################

# Function to get spectral feature columns based on platform
get_spectral_features <- function(data, platform) {
  
  # Common vegetation indices
  vi_features <- c("NDVI", "GNDVI", "NDRE", "GLI", "NGRDI", "SAVI", 
                   "EVI", "MSAVI", "SIPI", "MTCI", "CIgreen", "CIrededge", 
                   "ARVI", "VARI", "OSAVI", "TGI", "ExG", "RGRI")
  
  if (grepl("Satellite", platform)) {
    # Satellite features
    band_features <- c("Red", "Green", "Blue", "NIR", "Red_edge", "Deep_blue")
    
    # Get all feature columns
    feature_cols <- c()
    
    # Band features
    for (band in band_features) {
      cols <- grep(paste0("Satellite_", band, "_"), names(data), value = TRUE)
      feature_cols <- c(feature_cols, cols)
    }
    
    # Vegetation indices
    for (vi in vi_features) {
      cols <- grep(paste0("^", vi, "_Satellite_"), names(data), value = TRUE)
      feature_cols <- c(feature_cols, cols)
    }
    
  } else {
    # UAV features - CORRECTED: Use MS_ prefix instead of UAV_
    band_features <- c("Red", "Green", "Blue", "NIR", "Red_edge")
    
    # Get all feature columns
    feature_cols <- c()
    
    # Band features with MS_ prefix
    for (band in band_features) {
      cols <- grep(paste0("MS_", band, "_"), names(data), value = TRUE)
      feature_cols <- c(feature_cols, cols)
    }
    
    # Vegetation indices with _MS_ suffix
    for (vi in vi_features) {
      cols <- grep(paste0("^", vi, "_MS_"), names(data), value = TRUE)
      feature_cols <- c(feature_cols, cols)
    }
  }
  
  # Return unique columns that exist in the data
  return(unique(feature_cols[feature_cols %in% names(data)]))
}

# Function to prepare stage-specific data from temporal data
prepare_stage_specific_temporal <- function(temporal_data, growth_stages = NULL, 
                                            platform = "Satellite") {
  
  log_message(sprintf("Preparing data for stages: %s", 
                      paste(growth_stages, collapse = ", ")))
  
  # Filter by growth stage if specified
  if (!is.null(growth_stages)) {
    temporal_data <- temporal_data %>%
      filter(Growth_Stage %in% growth_stages)
  }
  
  # Get spectral features for this platform
  spectral_features <- get_spectral_features(temporal_data, platform)
  log_message(sprintf("  Found %d spectral features", length(spectral_features)))
  
  # Group by plot and aggregate features
  stage_data <- temporal_data %>%
    group_by(Plot_ID, Location, Line, Leaf_type, Spacing) %>%
    summarise(
      # Target variable - take first (should be same for all timepoints)
      peak_lai = first(peak_lai),
      yield = first(yield),
      
      # Number of observations
      n_observations = n(),
      
      # For each spectral feature, calculate statistics
      across(all_of(spectral_features), 
             list(mean = ~mean(.x, na.rm = TRUE),
                  max = ~max(.x, na.rm = TRUE),
                  min = ~min(.x, na.rm = TRUE),
                  sd = ~sd(.x, na.rm = TRUE),
                  cv = ~sd(.x, na.rm = TRUE)/mean(.x, na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      
      .groups = "drop"
    )
  
  # Handle infinite values and NAs
  numeric_cols <- names(stage_data)[sapply(stage_data, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("Plot_ID", "peak_lai", "yield"))
  
  stage_data[numeric_cols] <- lapply(stage_data[numeric_cols], function(x) {
    x[is.infinite(x)] <- NA
    if (sum(!is.na(x)) > 0) {
      x[is.na(x)] <- median(x, na.rm = TRUE)
    } else {
      x[is.na(x)] <- 0
    }
    return(x)
  })
  
  # Remove columns with zero variance
  zero_var_cols <- numeric_cols[sapply(stage_data[numeric_cols], var, na.rm = TRUE) == 0]
  if (length(zero_var_cols) > 0) {
    log_message(sprintf("  Removing %d zero-variance columns", length(zero_var_cols)))
    stage_data <- stage_data[, !names(stage_data) %in% zero_var_cols]
  }
  
  log_message(sprintf("  Final dataset: %d plots, %d features", 
                      nrow(stage_data), ncol(stage_data) - 6))
  
  return(stage_data)
}

################################################################################
# SECTION 3: PLOT-AWARE SAMPLING FUNCTIONS (Same as original)
################################################################################

#' Create plot-aware stratified split with proper validation
create_plot_aware_split <- function(data, train_prop = 0.70, val_prop = 0.15, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Get unique plots with stratification variables
  plot_info <- data %>%
    group_by(Plot_ID) %>%
    summarise(
      peak_lai = first(peak_lai),
      Spacing = first(Spacing),
      Leaf_type = first(Leaf_type),
      .groups = "drop"
    )
  
  # Create stratification bins
  n_unique_lai <- length(unique(plot_info$peak_lai))
  n_bins <- min(5, n_unique_lai)
  
  if (n_bins > 1) {
    plot_info <- plot_info %>%
      mutate(
        lai_bin = cut(peak_lai, breaks = n_bins, labels = FALSE, include.lowest = TRUE),
        strata = paste(Spacing, Leaf_type, lai_bin, sep = "_")
      )
  } else {
    plot_info <- plot_info %>%
      mutate(
        lai_bin = 1,
        strata = paste(Spacing, Leaf_type, sep = "_")
      )
  }
  
  # Stratified sampling
  test_plots <- plot_info %>%
    group_by(strata) %>%
    mutate(n_strata = n()) %>%
    ungroup() %>%
    group_by(strata) %>%
    sample_n(size = pmax(1, round(n() * (1 - train_prop - val_prop)))) %>%
    ungroup() %>%
    pull(Plot_ID)
  
  remaining_plots <- setdiff(plot_info$Plot_ID, test_plots)
  
  val_plots <- plot_info %>%
    filter(Plot_ID %in% remaining_plots) %>%
    group_by(strata) %>%
    mutate(n_strata = n()) %>%
    ungroup() %>%
    group_by(strata) %>%
    sample_n(size = pmax(1, round(n() * val_prop / (train_prop + val_prop)))) %>%
    ungroup() %>%
    pull(Plot_ID)
  
  train_plots <- setdiff(remaining_plots, val_plots)
  
  # Create splits
  train_data <- data %>% filter(Plot_ID %in% train_plots)
  val_data <- data %>% filter(Plot_ID %in% val_plots)
  test_data <- data %>% filter(Plot_ID %in% test_plots)
  
  # Verify no overlap
  if (length(intersect(train_plots, test_plots)) > 0 ||
      length(intersect(train_plots, val_plots)) > 0 ||
      length(intersect(val_plots, test_plots)) > 0) {
    stop("Plot overlap detected in train/val/test split!")
  }
  
  return(list(
    train = train_data,
    val = val_data,
    test = test_data,
    train_plots = train_plots,
    val_plots = val_plots,
    test_plots = test_plots
  ))
}

#' Plot-aware bootstrap sampling
bootstrap_plot_sample <- function(data, size = NULL) {
  
  if (is.null(size)) size <- nrow(data)
  
  # Get unique plots
  unique_plots <- unique(data$Plot_ID)
  n_plots <- length(unique_plots)
  
  # Sample plots with replacement
  sampled_plots <- sample(unique_plots, size = n_plots, replace = TRUE)
  
  # Build bootstrap sample
  boot_data <- data.frame()
  for (plot in sampled_plots) {
    plot_data <- data %>% filter(Plot_ID == plot)
    boot_data <- rbind(boot_data, plot_data)
  }
  
  return(boot_data)
}

################################################################################
# SECTION 4: MODEL TRAINING FUNCTIONS (Same as original)
################################################################################

#' Train a single model
train_single_model <- function(model_name, X_train, y_train, X_val, y_val, 
                               X_test, y_test, feature_names) {
  
  tryCatch({
    
    if (model_name == "lasso") {
      cv_model <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = 5)
      model <- glmnet(X_train, y_train, alpha = 1, lambda = cv_model$lambda.min)
      pred_test <- predict(model, X_test)[, 1]
      
      # Get non-zero coefficients for importance
      coefs <- coef(model)[-1, 1]  # Remove intercept
      importance <- abs(coefs)
      names(importance) <- feature_names
      
    } else if (model_name == "ridge") {
      cv_model <- cv.glmnet(X_train, y_train, alpha = 0, nfolds = 5)
      model <- glmnet(X_train, y_train, alpha = 0, lambda = cv_model$lambda.min)
      pred_test <- predict(model, X_test)[, 1]
      
      coefs <- coef(model)[-1, 1]
      importance <- abs(coefs)
      names(importance) <- feature_names
      
    } else if (model_name == "lm") {
      # Select top features for linear model
      cors <- abs(cor(X_train, y_train))
      top_features <- order(cors, decreasing = TRUE)[1:min(30, ncol(X_train))]
      
      train_df <- data.frame(y = y_train, X_train[, top_features])
      test_df <- data.frame(X_test[, top_features])
      
      model <- lm(y ~ ., data = train_df)
      pred_test <- predict(model, test_df)
      importance <- NULL
      
    } else if (model_name == "pls") {
      model <- plsr(y_train ~ X_train, validation = "CV", segments = 5)
      ncomp <- selectNcomp(model, method = "onesigma")
      pred_test <- predict(model, X_test, ncomp = ncomp)[, 1, 1]
      importance <- NULL
      
    } else if (model_name == "rf") {
      model <- randomForest(X_train, y_train, ntree = 500, mtry = sqrt(ncol(X_train)),
                            importance = TRUE)
      pred_test <- predict(model, X_test)
      importance <- model$importance[, "%IncMSE"]
      
    } else if (model_name == "xgboost") {
      dtrain <- xgb.DMatrix(X_train, label = y_train)
      dval <- xgb.DMatrix(X_val, label = y_val)
      
      params <- list(
        objective = "reg:squarederror",
        max_depth = 6,
        eta = 0.1,
        subsample = 0.8,
        colsample_bytree = 0.8
      )
      
      model <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = 100,
        watchlist = list(val = dval),
        early_stopping_rounds = 10,
        verbose = 0
      )
      
      pred_test <- predict(model, xgb.DMatrix(X_test))
      
      # Get feature importance
      imp_matrix <- xgb.importance(model = model)
      importance <- setNames(imp_matrix$Gain, imp_matrix$Feature)
      
    } else if (model_name == "svm") {
      model <- svm(X_train, y_train, kernel = "radial", gamma = 1/ncol(X_train))
      pred_test <- predict(model, X_test)
      importance <- NULL
    }
    
    # Calculate metrics
    metrics <- calculate_metrics_enhanced(y_test, pred_test)
    
    return(list(
      predictions = pred_test,
      metrics = metrics,
      importance = importance
    ))
    
  }, error = function(e) {
    warning(sprintf("Error training %s: %s", model_name, e$message))
    return(NULL)
  })
}

#' Calculate enhanced metrics
calculate_metrics_enhanced <- function(observed, predicted) {
  
  valid_idx <- !is.na(observed) & !is.na(predicted)
  obs <- observed[valid_idx]
  pred <- predicted[valid_idx]
  
  if (length(obs) < 5) {
    return(list(R2 = NA, RMSE = NA, MAE = NA, Bias = NA, n = length(obs)))
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

#' Aggregate results across repetitions
aggregate_repetition_results <- function(rep_results, models, stage_name, 
                                         dataset_name, n_reps) {
  
  if (length(rep_results) == 0) {
    log_message("  No results to aggregate", "WARNING")
    return(NULL)
  }
  
  aggregated <- list()
  
  for (model_name in models) {
    # Extract metrics for this model across all repetitions
    model_metrics <- list()
    
    for (rep in 1:length(rep_results)) {
      if (!is.null(rep_results[[rep]]$metrics[[model_name]])) {
        model_metrics[[rep]] <- rep_results[[rep]]$metrics[[model_name]]
      }
    }
    
    if (length(model_metrics) > 0) {
      # Calculate statistics
      metrics_df <- do.call(rbind, lapply(model_metrics, as.data.frame))
      
      # Calculate mean, SD, and CI for each metric
      summary_stats <- list()
      
      for (metric in c("R2", "RMSE", "MAE", "Bias", "MAPE")) {
        values <- metrics_df[[metric]]
        values <- values[!is.na(values)]
        
        if (length(values) > 0) {
          summary_stats[[metric]] <- list(
            mean = mean(values),
            sd = ifelse(length(values) > 1, sd(values), 0),
            ci_lower = quantile(values, (1 - CONFIDENCE_LEVEL) / 2),
            ci_upper = quantile(values, 1 - (1 - CONFIDENCE_LEVEL) / 2),
            n = length(values)
          )
        }
      }
      
      aggregated[[model_name]] <- list(
        metrics = summary_stats,
        n_successful_reps = length(model_metrics),
        dataset = dataset_name,
        stage = stage_name
      )
    }
  }
  
  # Aggregate predictions
  all_predictions <- data.frame()
  for (rep in 1:length(rep_results)) {
    for (model_name in names(rep_results[[rep]]$predictions)) {
      all_predictions <- rbind(all_predictions, 
                               rep_results[[rep]]$predictions[[model_name]])
    }
  }
  
  aggregated$all_predictions <- all_predictions
  aggregated$n_repetitions <- n_reps
  
  return(aggregated)
}

# Setup parallel processing
if (USE_PARALLEL) {
  n_cores <- min(detectCores() - 1, 32)  # Limit to 32 cores on cluster
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export all necessary functions and libraries to cluster workers
  clusterExport(cl, c("create_plot_aware_split", "bootstrap_plot_sample", 
                      "train_single_model", "calculate_metrics_enhanced",
                      "log_message", "CONFIDENCE_LEVEL"))
  
  # Load libraries on workers
  clusterEvalQ(cl, {
    library(tidyverse)
    library(glmnet)
    library(randomForest)
    library(xgboost)
    library(e1071)
    library(pls)
  })
} else {
  registerDoSEQ()
}

################################################################################
# SECTION 5: ENHANCED MODEL TRAINING WITH CONFIDENCE INTERVALS
################################################################################

#' Train models with bootstrap confidence intervals
train_models_with_ci <- function(data, stage_name, dataset_name, 
                                 n_reps = N_REPETITIONS, 
                                 n_bootstrap = N_BOOTSTRAP) {
  
  log_message(sprintf("Training models for %s - %s with %d repetitions", 
                      dataset_name, stage_name, n_reps))
  
  # Check data size
  if (nrow(data) < 30) {
    log_message(sprintf("  ERROR: Insufficient data for modeling (%d rows)", nrow(data)), "ERROR")
    return(NULL)
  }
  
  # Check for valid numeric features
  numeric_cols <- sapply(data, is.numeric)
  exclude_cols <- c("Plot_ID", "peak_lai", "yield", "Spacing", "n_observations")
  feature_cols <- names(data)[numeric_cols & !names(data) %in% exclude_cols]
  
  if (length(feature_cols) < 5) {
    log_message(sprintf("  ERROR: Insufficient numeric features (%d features)", 
                        length(feature_cols)), "ERROR")
    return(NULL)
  }
  
  log_message(sprintf("  Data shape: %d rows x %d features", nrow(data), length(feature_cols)))
  
  # Define models
  models <- c("lasso", "ridge", "lm", "pls", "rf", "xgboost", "svm")
  
  # Storage for all repetitions
  all_results <- list()
  all_predictions <- list()
  all_feature_importance <- list()
  
  # Progress tracking
  pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
  
  # Parallel execution of repetitions
  rep_results <- foreach(rep = 1:n_reps, 
                         .packages = c("tidyverse", "glmnet", "randomForest", 
                                       "xgboost", "e1071", "pls", "caret"),
                         .errorhandling = 'pass') %dopar% {
                           
                           tryCatch({
                             # Create plot-aware split for this repetition
                             splits <- create_plot_aware_split(data, seed = rep * 1000)
                             
                             # Check split sizes
                             if (nrow(splits$train) < 20) {
                               return(list(error = sprintf("Training set too small: %d rows", 
                                                           nrow(splits$train))))
                             }
                             if (nrow(splits$val) < 5) {
                               return(list(error = sprintf("Validation set too small: %d rows", 
                                                           nrow(splits$val))))
                             }
                             if (nrow(splits$test) < 5) {
                               return(list(error = sprintf("Test set too small: %d rows", 
                                                           nrow(splits$test))))
                             }
                             
                             # Prepare features
                             X_train <- as.matrix(splits$train[, feature_cols])
                             y_train <- splits$train$peak_lai
                             X_val <- as.matrix(splits$val[, feature_cols])
                             y_val <- splits$val$peak_lai
                             X_test <- as.matrix(splits$test[, feature_cols])
                             y_test <- splits$test$peak_lai
                             
                             # Results for this repetition
                             rep_metrics <- list()
                             rep_predictions <- list()
                             rep_importance <- list()
                             
                             # Train each model
                             for (model_name in models) {
                               
                               model_results <- train_single_model(
                                 model_name, X_train, y_train, X_val, y_val, X_test, y_test,
                                 feature_names = feature_cols
                               )
                               
                               if (!is.null(model_results)) {
                                 rep_metrics[[model_name]] <- model_results$metrics
                                 rep_predictions[[model_name]] <- data.frame(
                                   Plot_ID = splits$test$Plot_ID,
                                   Observed = y_test,
                                   Predicted = model_results$predictions,
                                   Model = model_name,
                                   Repetition = rep
                                 )
                                 
                                 # Store feature importance if available
                                 if (!is.null(model_results$importance)) {
                                   rep_importance[[model_name]] <- model_results$importance
                                 }
                               }
                             }
                             
                             list(
                               metrics = rep_metrics,
                               predictions = rep_predictions,
                               importance = rep_importance,
                               n_train = nrow(splits$train),
                               n_val = nrow(splits$val),
                               n_test = nrow(splits$test)
                             )
                             
                           }, error = function(e) {
                             error_msg <- paste("Error in repetition", rep, ":", e$message)
                             warning(error_msg)
                             return(list(error = e$message))
                           })
                         }
  
  close(pb)
  
  # Check for errors in results
  errors <- sapply(rep_results, function(x) !is.null(x$error))
  if (any(errors)) {
    log_message("  Errors encountered in some repetitions:", "WARNING")
    error_msgs <- unique(unlist(lapply(rep_results[errors], function(x) x$error)))
    for (msg in error_msgs) {
      log_message(sprintf("    - %s", msg), "WARNING")
    }
  }
  
  # Remove failed repetitions
  rep_results <- rep_results[!sapply(rep_results, is.null) & !errors]
  n_successful <- length(rep_results)
  
  log_message(sprintf("  Successfully completed %d/%d repetitions", n_successful, n_reps))
  
  if (n_successful == 0) {
    log_message("  ERROR: All repetitions failed", "ERROR")
    return(NULL)
  }
  
  # Aggregate results across repetitions
  aggregated_results <- aggregate_repetition_results(
    rep_results, models, stage_name, dataset_name, n_successful
  )
  
  # Bootstrap feature importance for models that support it
  if (n_bootstrap > 0 && length(feature_cols) > 0) {
    log_message("  Calculating bootstrap feature importance...")
    aggregated_results$bootstrap_importance <- calculate_bootstrap_importance(
      data, feature_cols, models = c("rf", "xgboost", "lasso", "ridge"),
      n_bootstrap = n_bootstrap
    )
  }
  
  return(aggregated_results)
}

#' Calculate bootstrap feature importance
calculate_bootstrap_importance <- function(data, feature_cols, models, n_bootstrap) {
  
  importance_results <- list()
  
  for (model_name in models) {
    if (model_name %in% c("rf", "xgboost", "lasso", "ridge")) {
      
      boot_importance <- foreach(b = 1:n_bootstrap, 
                                 .packages = c("randomForest", "xgboost", "glmnet", "tidyverse"),
                                 .combine = 'rbind') %dopar% {
                                   
                                   # Bootstrap sample
                                   boot_data <- bootstrap_plot_sample(data)
                                   
                                   X_boot <- as.matrix(boot_data[, feature_cols])
                                   y_boot <- boot_data$peak_lai
                                   
                                   if (model_name == "rf") {
                                     model <- randomForest(X_boot, y_boot, ntree = 200, importance = TRUE)
                                     importance <- model$importance[, "%IncMSE"]
                                     
                                   } else if (model_name == "xgboost") {
                                     dtrain <- xgb.DMatrix(X_boot, label = y_boot)
                                     model <- xgb.train(
                                       params = list(objective = "reg:squarederror", max_depth = 6, eta = 0.1),
                                       data = dtrain,
                                       nrounds = 50,
                                       verbose = 0
                                     )
                                     imp_matrix <- xgb.importance(model = model)
                                     importance <- setNames(imp_matrix$Gain, imp_matrix$Feature)
                                     
                                     # Match to all features
                                     full_importance <- rep(0, length(feature_cols))
                                     names(full_importance) <- feature_cols
                                     full_importance[names(importance)] <- importance
                                     importance <- full_importance
                                     
                                   } else if (model_name %in% c("lasso", "ridge")) {
                                     alpha <- ifelse(model_name == "lasso", 1, 0)
                                     cv_model <- cv.glmnet(X_boot, y_boot, alpha = alpha, nfolds = 5)
                                     model <- glmnet(X_boot, y_boot, alpha = alpha, lambda = cv_model$lambda.min)
                                     coefs <- coef(model)[-1, 1]
                                     importance <- abs(coefs)
                                     names(importance) <- feature_cols
                                   }
                                   
                                   importance
                                 }
      
      # Aggregate importance across bootstraps
      if (!is.null(boot_importance)) {
        if (is.matrix(boot_importance)) {
          imp_summary <- data.frame(
            feature = colnames(boot_importance),
            importance_mean = colMeans(boot_importance, na.rm = TRUE),
            importance_sd = apply(boot_importance, 2, sd, na.rm = TRUE),
            importance_ci_lower = apply(boot_importance, 2, quantile, 
                                        probs = (1 - CONFIDENCE_LEVEL) / 2, na.rm = TRUE),
            importance_ci_upper = apply(boot_importance, 2, quantile, 
                                        probs = 1 - (1 - CONFIDENCE_LEVEL) / 2, na.rm = TRUE)
          ) %>%
            arrange(desc(importance_mean))
        } else {
          # Handle vector case
          imp_summary <- data.frame(
            feature = names(boot_importance),
            importance_mean = boot_importance,
            importance_sd = 0,
            importance_ci_lower = boot_importance,
            importance_ci_upper = boot_importance
          ) %>%
            arrange(desc(importance_mean))
        }
        
        importance_results[[model_name]] <- imp_summary
      }
    }
  }
  
  return(importance_results)
}

################################################################################
# SECTION 6: STATISTICAL TESTING FUNCTIONS (Same as original)
################################################################################

#' Perform pairwise statistical tests between models
compare_models_statistically <- function(results_list, metric = "R2") {
  
  model_names <- names(results_list)
  n_models <- length(model_names)
  
  # Initialize p-value matrix
  p_matrix <- matrix(1, nrow = n_models, ncol = n_models)
  rownames(p_matrix) <- colnames(p_matrix) <- model_names
  
  # Extract metric values for each model
  metric_values <- list()
  for (model in model_names) {
    if (!is.null(results_list[[model]]$all_predictions)) {
      predictions <- results_list[[model]]$all_predictions
      
      # Calculate metric for each repetition
      rep_metrics <- predictions %>%
        group_by(Repetition) %>%
        summarise(
          value = calculate_metrics_enhanced(Observed, Predicted)[[metric]],
          .groups = "drop"
        ) %>%
        pull(value)
      
      metric_values[[model]] <- rep_metrics
    }
  }
  
  # Pairwise comparisons
  for (i in 1:(n_models-1)) {
    for (j in (i+1):n_models) {
      if (length(metric_values[[i]]) > 0 && length(metric_values[[j]]) > 0) {
        # Wilcoxon signed-rank test
        test_result <- wilcox.test(metric_values[[i]], metric_values[[j]], 
                                   paired = FALSE, exact = FALSE)
        p_matrix[i, j] <- p_matrix[j, i] <- test_result$p.value
      }
    }
  }
  
  return(p_matrix)
}

################################################################################
# SECTION 7: STAGE-SPECIFIC ANALYSIS WITH TEMPORAL DATA
################################################################################

log_message("\n=== STAGE-SPECIFIC MODEL ANALYSIS WITH TEMPORAL DATA ===")

# Get unique growth stages from temporal data
all_stages <- unique(unlist(lapply(temporal_datasets, function(x) unique(x$Growth_Stage))))
log_message(sprintf("Available growth stages: %s", paste(all_stages, collapse = ", ")))

# Define stage combinations
stage_combinations <- list(
  "Vegetative" = c("Vegetative"),
  "Reproductive" = c("Reproductive"),
  "Maturity" = c("Maturity"),
  "Veg_Rep" = c("Vegetative", "Reproductive"),
  "Veg_Mat" = c("Vegetative", "Maturity"),
  "Rep_Mat" = c("Reproductive", "Maturity"),
  "All_Stages" = c("Vegetative", "Reproductive", "Maturity")
)

# Filter combinations based on available stages
stage_combinations <- stage_combinations[sapply(stage_combinations, function(x) all(x %in% all_stages))]
log_message(sprintf("Using %d stage combinations", length(stage_combinations)))

stage_results <- list()
all_bootstrap_results <- list()

for (dataset_name in names(temporal_datasets)) {
  log_message(sprintf("\nProcessing %s", dataset_name))
  temporal_data <- temporal_datasets[[dataset_name]]
  
  # Determine platform
  platform <- ifelse(grepl("Satellite", dataset_name), "Satellite", "UAV")
  
  for (stage_name in names(stage_combinations)) {
    log_message(sprintf("\n  Analyzing %s stage combination", stage_name))
    
    # Get the stages to include
    stages_to_include <- stage_combinations[[stage_name]]
    
    # Prepare data for this stage combination
    stage_data <- prepare_stage_specific_temporal(
      temporal_data, 
      growth_stages = stages_to_include,
      platform = platform
    )
    
    if (!is.null(stage_data) && nrow(stage_data) > 50) {
      # Train models with confidence intervals
      model_results <- train_models_with_ci(
        stage_data, 
        stage_name, 
        dataset_name,
        n_reps = N_REPETITIONS,
        n_bootstrap = N_BOOTSTRAP
      )
      
      if (length(model_results) > 0) {
        stage_results[[dataset_name]][[stage_name]] <- model_results
        
        # Store bootstrap results
        if (!is.null(model_results$bootstrap_importance)) {
          all_bootstrap_results[[dataset_name]][[stage_name]] <- 
            model_results$bootstrap_importance
        }
      }
    } else {
      log_message(sprintf("    Insufficient data for %s", stage_name), "WARNING")
    }
  }
}

# Save enhanced results
saveRDS(stage_results, file.path(output_dir, "stage_specific_models", 
                                 "all_results_temporal_with_ci.rds"))
saveRDS(all_bootstrap_results, file.path(output_dir, "bootstrap_results",
                                         "feature_importance_temporal_bootstrap.rds"))

################################################################################
# SECTION 8: COMPILE ENHANCED PERFORMANCE SUMMARY
################################################################################

log_message("\n=== Compiling Enhanced Performance Summary ===")

performance_summary_ci <- data.frame()

for (dataset_name in names(stage_results)) {
  for (stage_name in names(stage_results[[dataset_name]])) {
    stage_result <- stage_results[[dataset_name]][[stage_name]]
    
    for (model_name in names(stage_result)) {
      if (model_name %in% c("all_predictions", "n_repetitions", "bootstrap_importance")) {
        next
      }
      
      model_result <- stage_result[[model_name]]
      
      if (!is.null(model_result$metrics)) {
        # Extract metrics with CI
        metrics <- model_result$metrics
        
        performance_summary_ci <- rbind(performance_summary_ci, data.frame(
          Dataset = dataset_name,
          Stage = stage_name,
          Model = model_name,
          R2_mean = metrics$R2$mean,
          R2_sd = metrics$R2$sd,
          R2_ci_lower = metrics$R2$ci_lower,
          R2_ci_upper = metrics$R2$ci_upper,
          RMSE_mean = metrics$RMSE$mean,
          RMSE_ci_lower = metrics$RMSE$ci_lower,
          RMSE_ci_upper = metrics$RMSE$ci_upper,
          MAE_mean = metrics$MAE$mean,
          MAE_ci_lower = metrics$MAE$ci_lower,
          MAE_ci_upper = metrics$MAE$ci_upper,
          MAPE_mean = if (!is.null(metrics$MAPE)) metrics$MAPE$mean else NA,
          Bias_mean = metrics$Bias$mean,
          N_repetitions = model_result$n_successful_reps,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Add formatted columns
if (nrow(performance_summary_ci) > 0) {
  performance_summary_ci <- performance_summary_ci %>%
    mutate(
      R2_formatted = sprintf("%.3f [%.3f, %.3f]", R2_mean, R2_ci_lower, R2_ci_upper),
      RMSE_formatted = sprintf("%.3f [%.3f, %.3f]", RMSE_mean, RMSE_ci_lower, RMSE_ci_upper),
      MAE_formatted = sprintf("%.3f [%.3f, %.3f]", MAE_mean, MAE_ci_lower, MAE_ci_upper)
    )
  
  write_csv(performance_summary_ci, 
            file.path(output_dir, "confidence_intervals", 
                      "stage_specific_performance_temporal_with_ci.csv"))
  
  log_message(sprintf("Performance summary created with %d entries", nrow(performance_summary_ci)))
} else {
  log_message("Warning: No performance results to summarize", "WARNING")
}

################################################################################
# SECTION 9: STATISTICAL COMPARISONS
################################################################################

log_message("\n=== Performing Statistical Comparisons ===")

# Compare models within each stage and dataset
statistical_comparisons <- list()

for (dataset_name in names(stage_results)) {
  for (stage_name in names(stage_results[[dataset_name]])) {
    
    log_message(sprintf("  Statistical tests for %s - %s", dataset_name, stage_name))
    
    # Get all model results for this stage
    stage_models <- stage_results[[dataset_name]][[stage_name]]
    
    # Remove non-model entries
    model_names <- setdiff(names(stage_models), 
                           c("all_predictions", "n_repetitions", "bootstrap_importance"))
    
    if (length(model_names) > 1) {
      # Extract model results
      model_list <- list()
      for (model in model_names) {
        model_list[[model]] <- stage_models[[model]]
      }
      
      # Perform pairwise comparisons
      p_matrix <- compare_models_statistically(model_list, metric = "R2")
      
      statistical_comparisons[[paste(dataset_name, stage_name, sep = "_")]] <- p_matrix
      
      # Save p-value matrix
      write.csv(p_matrix, 
                file.path(output_dir, "statistical_tests",
                          sprintf("pvalues_temporal_%s_%s.csv", dataset_name, stage_name)))
    }
  }
}

################################################################################
# SECTION 10: ENHANCED VISUALIZATIONS
################################################################################

log_message("\n=== Creating Enhanced Visualizations ===")

# 1. Performance comparison with error bars
create_performance_comparison_plot <- function(performance_data, output_file) {
  
  p <- ggplot(performance_data, aes(x = Stage, y = R2_mean, color = Model)) +
    geom_point(position = position_dodge(0.5), size = 3) +
    geom_errorbar(aes(ymin = R2_ci_lower, ymax = R2_ci_upper), 
                  width = 0.2, position = position_dodge(0.5)) +
    facet_wrap(~ Dataset, scales = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") +
    labs(title = "Model Performance with 95% Confidence Intervals (Temporal Data)",
         x = "Growth Stage Combination", 
         y = "R² (mean with 95% CI)") +
    scale_color_brewer(palette = "Set2")
  
  ggsave(output_file, plot = p, width = 14, height = 8)
  return(p)
}

# Create performance plot
if (nrow(performance_summary_ci) > 0) {
  create_performance_comparison_plot(
    performance_summary_ci,
    file.path(output_dir, "confidence_intervals", "performance_temporal_with_ci.png")
  )
}

################################################################################
# SECTION 11: TEMPORAL ACCUMULATION ANALYSIS
################################################################################

log_message("\n=== Temporal Accumulation Analysis ===")

# Test progressive prediction
temporal_accumulation_results_ci <- data.frame()

for (dataset_name in names(temporal_datasets)) {
  log_message(sprintf("\nTemporal accumulation for %s", dataset_name))
  
  # Define accumulation sequence
  accumulation_stages <- list(
    "Vegetative_Only" = c("Vegetative"),
    "Veg_Plus_Reproductive" = c("Vegetative", "Reproductive"),
    "All_Stages" = c("Vegetative", "Reproductive", "Maturity")
  )
  
  for (acc_name in names(accumulation_stages)) {
    stages <- accumulation_stages[[acc_name]]
    
    # Get results for this combination
    stage_key <- paste(stages, collapse = "_")
    if (stage_key == "Vegetative") stage_key <- "Vegetative"
    else if (stage_key == "Vegetative_Reproductive") stage_key <- "Veg_Rep"
    else if (stage_key == "Vegetative_Reproductive_Maturity") stage_key <- "All_Stages"
    
    if (!is.null(stage_results[[dataset_name]][[stage_key]])) {
      for (model_name in names(stage_results[[dataset_name]][[stage_key]])) {
        if (model_name %in% c("all_predictions", "n_repetitions", "bootstrap_importance")) {
          next
        }
        
        model_result <- stage_results[[dataset_name]][[stage_key]][[model_name]]
        
        if (!is.null(model_result$metrics)) {
          temporal_accumulation_results_ci <- rbind(temporal_accumulation_results_ci, 
                                                    data.frame(
                                                      Dataset = dataset_name,
                                                      Accumulation = acc_name,
                                                      Model = model_name,
                                                      R2_mean = model_result$metrics$R2$mean,
                                                      R2_ci_lower = model_result$metrics$R2$ci_lower,
                                                      R2_ci_upper = model_result$metrics$R2$ci_upper,
                                                      RMSE_mean = model_result$metrics$RMSE$mean,
                                                      RMSE_ci_lower = model_result$metrics$RMSE$ci_lower,
                                                      RMSE_ci_upper = model_result$metrics$RMSE$ci_upper,
                                                      stringsAsFactors = FALSE
                                                    )
          )
        }
      }
    }
  }
}

write_csv(temporal_accumulation_results_ci, 
          file.path(output_dir, "temporal_accumulation", 
                    "accumulation_results_temporal_with_ci.csv"))

################################################################################
# SECTION 12: ENHANCED SUMMARY REPORT
################################################################################

log_message("\n=== Creating Enhanced Summary Report ===")

# Create comprehensive summary
sink(file.path(output_dir, "phase4_growth_stage_peak_lai_temporal_summary.txt"))

cat("GROWTH STAGE ANALYSIS FOR PEAK LAI - TEMPORAL DATA VERSION WITH CONFIDENCE INTERVALS\n")
cat("===================================================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("ANALYSIS CONFIGURATION:\n")
cat(sprintf("- Repeated holdout validation: %d repetitions\n", N_REPETITIONS))
cat(sprintf("- Bootstrap iterations for feature importance: %d\n", N_BOOTSTRAP))
cat(sprintf("- Confidence level: %.0f%%\n", CONFIDENCE_LEVEL * 100))
cat(sprintf("- Train/Val/Test split: 70%%/15%%/15%%\n"))
cat("- Plot-aware sampling to prevent data leakage\n")
cat("- Using raw temporal data with dynamic feature aggregation\n\n")

cat("KEY FINDINGS:\n\n")

if (nrow(performance_summary_ci) > 0) {
  # Best performing model per stage with CI
  cat("1. BEST MODELS BY STAGE (with 95% CI):\n")
  best_models <- performance_summary_ci %>%
    group_by(Stage) %>%
    slice_max(R2_mean, n = 1) %>%
    select(Stage, Model, R2_formatted, RMSE_formatted)
  
  for (i in 1:nrow(best_models)) {
    cat(sprintf("   %s: %s model, R² = %s, RMSE = %s\n",
                best_models$Stage[i], 
                best_models$Model[i],
                best_models$R2_formatted[i],
                best_models$RMSE_formatted[i]))
  }
  
  cat("\n2. STAGE COMPARISON:\n")
  # Average performance by stage
  stage_avg <- performance_summary_ci %>%
    group_by(Stage) %>%
    summarise(
      Mean_R2 = mean(R2_mean, na.rm = TRUE),
      SD_R2 = sd(R2_mean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_R2))
  
  cat("   Average R² by stage:\n")
  for (i in 1:nrow(stage_avg)) {
    cat(sprintf("   %s: %.3f ± %.3f\n",
                stage_avg$Stage[i],
                stage_avg$Mean_R2[i],
                stage_avg$SD_R2[i]))
  }
}

cat("\n3. DATA SUMMARY:\n")
for (dataset in names(temporal_datasets)) {
  n_obs <- nrow(temporal_datasets[[dataset]])
  n_plots <- length(unique(temporal_datasets[[dataset]]$Plot_ID))
  stages <- unique(temporal_datasets[[dataset]]$Growth_Stage)
  cat(sprintf("   %s: %d observations, %d plots, stages: %s\n",
              dataset, n_obs, n_plots, paste(stages, collapse = ", ")))
}

sink()

################################################################################
# SECTION 13: CLEANUP
################################################################################

log_message("\n=== Phase 4 Growth Stage Analysis Complete (Temporal Version) ===")
log_message(sprintf("Results saved to: %s", output_dir))
log_message("Analysis used raw temporal data with dynamic feature aggregation")

# Clean up parallel cluster if used
if (USE_PARALLEL && exists("cl")) {
  stopCluster(cl)
}

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), 
           file.path(output_dir, "session_info_temporal.txt"))