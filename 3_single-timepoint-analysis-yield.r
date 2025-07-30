################################################################################
# Phase 2 Yield Prediction with Enhanced Statistical Rigor - CLUSTER VERSION
# 
# Purpose: Predict yield from spectral data with comprehensive statistical
#          uncertainty quantification matching Phase 4 methodology
#          
# Key Enhancements:
# - Plot-aware bootstrap resampling (100 iterations)
# - Repeated holdout validation (50 repetitions)
# - 95% confidence intervals for all metrics
# - Statistical significance testing between models
# - Uncertainty in feature importance
# - Parallel processing disabled for cluster stability
# - Both individual timepoints and temporal aggregation
# - INCREMENTAL SAVING AND ERROR RECOVERY
#
# Models: lasso, ridge, lm, pls, rf, xgboost, svm
#
# Author: Soybean Yield Prediction Pipeline - Enhanced
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
USE_PARALLEL <- FALSE  # Set to FALSE for cluster stability
N_BOOTSTRAP <- ifelse(DEBUG_MODE, 10, 100)  # Number of bootstrap iterations
N_REPETITIONS <- ifelse(DEBUG_MODE, 5, 50)  # Number of repeated holdout splits
CONFIDENCE_LEVEL <- 0.95  # For confidence intervals

# Setup parallel processing (disabled for cluster)
if (USE_PARALLEL) {
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
} else {
  registerDoSEQ()
}

# Check if on cluster
on_cluster <- Sys.getenv("SLURM_JOB_ID") != ""
if (on_cluster) {
  cat(sprintf("Running on cluster - Job ID: %s\n", Sys.getenv("SLURM_JOB_ID")))
}

################################################################################
# SECTION 0: SETUP AND PATHS - MODIFIED FOR CLUSTER STRUCTURE
################################################################################

# Set paths - MODIFIED for cluster directory structure
current_dir <- getwd()
project_root <- ifelse(basename(current_dir) == "scripts", dirname(current_dir), current_dir)

# Direct paths without phase subdirectories
data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "soybean_analysis_outputs", "phase2_yield_enhanced")

# Create directories
subdirs <- c("models", "predictions", "bootstrap_results", "confidence_intervals", 
             "statistical_tests", "visualizations", "diagnostics", "temporal_analysis",
             "individual_timepoints", "checkpoints")
for (subdir in subdirs) {
  dir_path <- file.path(output_dir, subdir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

# Logging function
log_file <- file.path(output_dir, "yield_analysis_enhanced_log.txt")
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- sprintf("[%s] %s: %s", timestamp, level, msg)
  cat(formatted_msg, "\n")
  cat(formatted_msg, "\n", file = log_file, append = TRUE)
}

log_message("Starting Phase 2: Yield Prediction with Enhanced Statistical Rigor - CLUSTER VERSION")
log_message(sprintf("Configuration: %d bootstrap iterations, %d repeated holdouts", 
                    N_BOOTSTRAP, N_REPETITIONS))

################################################################################
# SECTION 1: DATA LOADING - MODIFIED TO READ FROM DATA FOLDER
################################################################################

log_message("Loading datasets from data folder...")

# Define datasets - using full spectral timepoint files
datasets <- list(
  SF_UAV = file.path(data_dir, "SF_UAV_with_DAP_cleaned.csv"),
  EF_UAV = file.path(data_dir, "EF_UAV_with_DAP_cleaned.csv"),
  SF_Satellite = file.path(data_dir, "SF_Satellite_with_DAP_cleaned.csv"),
  EF_Satellite = file.path(data_dir, "EF_Satellite_with_DAP_cleaned.csv")
)

# Check if files exist
for (dataset_name in names(datasets)) {
  if (!file.exists(datasets[[dataset_name]])) {
    log_message(sprintf("ERROR: %s not found at %s", dataset_name, datasets[[dataset_name]]), "ERROR")
    stop(sprintf("Required file %s not found. Please ensure all DAP-cleaned files are in the data folder.", dataset_name))
  }
}

# Load data
loaded_data <- list()
for (dataset_name in names(datasets)) {
  file_path <- datasets[[dataset_name]]
  loaded_data[[dataset_name]] <- read_csv(file_path, show_col_types = FALSE)
  log_message(sprintf("Loaded %s: %d rows", dataset_name, nrow(loaded_data[[dataset_name]])))
}

################################################################################
# SECTION 2: DATA PREPARATION FUNCTIONS
################################################################################

#' Prepare yield data for individual timepoints
prepare_yield_timepoint_data <- function(data, timepoint) {
  
  # Filter to specific timepoint
  tp_data <- data %>%
    filter(Timepoint == timepoint) %>%
    filter(!grepl("Secondary", Line, ignore.case = TRUE))
  
  if (nrow(tp_data) == 0) {
    return(NULL)
  }
  
  # Select spectral features
  feature_patterns <- "Mean|Median|P10|P25|P75|P90"
  spectral_patterns <- "MS_|NDVI|GNDVI|NDRE|NGRDI|SAVI|MSAVI|MTCI|CIgreen|CIrededge|OSAVI|RGRI|Satellite_|EVI|ARVI|VARI|TGI|ExG|GLI|SIPI"
  
  feature_cols <- names(tp_data)[grepl(feature_patterns, names(tp_data))]
  spectral_cols <- names(tp_data)[grepl(spectral_patterns, names(tp_data))]
  feature_cols <- unique(c(feature_cols, spectral_cols))
  
  # Exclude metadata
  exclude_patterns <- c("Plot_ID", "Location", "Timepoint", "Folder", "Line", 
                        "Leaf_type", "Spacing", "yield", "peak_lai", "LAI_TP",
                        "Acquisition_Date", "ProcessedFile", "Platform", 
                        "Pass", "Range", "ID", "Location.x", "Location.y", 
                        "DAP", "Growth_Stage")
  
  for (pattern in exclude_patterns) {
    feature_cols <- feature_cols[!grepl(pattern, feature_cols, ignore.case = TRUE)]
  }
  
  # Keep only numeric features
  numeric_features <- feature_cols[sapply(tp_data[feature_cols], is.numeric)]
  
  # Create modeling dataset
  model_data <- tp_data %>%
    select(Plot_ID, all_of(numeric_features), yield, Spacing, Leaf_type, DAP) %>%
    filter(!is.na(yield))
  
  # Handle duplicate plots
  duplicate_plots <- model_data %>%
    count(Plot_ID) %>%
    filter(n > 1)
  
  if (nrow(duplicate_plots) > 0) {
    model_data <- model_data %>%
      group_by(Plot_ID) %>%
      slice_head(n = 1) %>%
      ungroup()
  }
  
  # Remove columns with >50% NA values
  model_data <- model_data %>%
    select_if(~sum(is.na(.)) < nrow(model_data) * 0.5)
  
  # Impute remaining NAs with median
  numeric_cols <- names(model_data)[sapply(model_data, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("Plot_ID", "yield", "DAP"))
  
  for (col in numeric_cols) {
    if (any(is.na(model_data[[col]]))) {
      model_data[[col]][is.na(model_data[[col]])] <- median(model_data[[col]], na.rm = TRUE)
    }
  }
  
  return(model_data)
}

#' Prepare yield data using temporal aggregation
prepare_yield_temporal_aggregation <- function(data) {
  
  log_message("Preparing yield data using temporal aggregation")
  
  # Filter out Secondary lines
  data <- data %>%
    filter(!grepl("Secondary", Line, ignore.case = TRUE))
  
  # Get spectral feature columns
  feature_patterns <- "Mean|Median|P10|P25|P75|P90"
  spectral_patterns <- "MS_|NDVI|GNDVI|NDRE|NGRDI|SAVI|MSAVI|MTCI|CIgreen|CIrededge|OSAVI|RGRI|Satellite_|EVI|ARVI|VARI|TGI|ExG|GLI|SIPI"
  
  feature_cols <- names(data)[grepl(feature_patterns, names(data))]
  spectral_cols <- names(data)[grepl(spectral_patterns, names(data))]
  feature_cols <- unique(c(feature_cols, spectral_cols))
  
  # Exclude metadata
  exclude_patterns <- c("Plot_ID", "Location", "Timepoint", "Folder", "Line", 
                        "Leaf_type", "Spacing", "yield", "peak_lai",
                        "Acquisition_Date", "ProcessedFile", "Platform", 
                        "Pass", "Range", "ID", "Location.x", "Location.y", 
                        "DAP", "Growth_Stage")
  
  for (pattern in exclude_patterns) {
    feature_cols <- feature_cols[!grepl(pattern, feature_cols, ignore.case = TRUE)]
  }
  
  # Keep only numeric features
  numeric_features <- feature_cols[sapply(data[feature_cols], is.numeric)]
  
  # Get unique yield values for each plot
  plot_yield <- data %>%
    select(Plot_ID, Spacing, Leaf_type, yield) %>%
    distinct() %>%
    filter(!is.na(yield))
  
  # Aggregate spectral features across all timepoints
  aggregated_features <- data %>%
    group_by(Plot_ID) %>%
    summarise(
      across(all_of(numeric_features),
             list(
               temporal_mean = ~mean(.x, na.rm = TRUE),
               temporal_max = ~max(.x, na.rm = TRUE),
               temporal_min = ~min(.x, na.rm = TRUE),
               temporal_sd = ~sd(.x, na.rm = TRUE),
               temporal_range = ~(max(.x, na.rm = TRUE) - min(.x, na.rm = TRUE)),
               temporal_cv = ~(sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE))
             ),
             .names = "{.col}_{.fn}"
      ),
      n_timepoints = n(),
      mean_dap = mean(DAP, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Join with yield and metadata
  aggregated_data <- aggregated_features %>%
    inner_join(plot_yield, by = "Plot_ID")
  
  # Remove any columns with all NA values
  aggregated_data <- aggregated_data %>%
    select_if(~!all(is.na(.)))
  
  # Impute remaining NAs with median
  numeric_cols <- names(aggregated_data)[sapply(aggregated_data, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("Plot_ID", "yield"))
  
  for (col in numeric_cols) {
    if (any(is.na(aggregated_data[[col]]))) {
      aggregated_data[[col]][is.na(aggregated_data[[col]])] <- median(aggregated_data[[col]], na.rm = TRUE)
    }
  }
  
  log_message(sprintf("Created %d temporal features for %d plots", 
                      ncol(aggregated_features) - 3, nrow(aggregated_data)))
  
  return(aggregated_data)
}

################################################################################
# SECTION 3: PLOT-AWARE SAMPLING FUNCTIONS
################################################################################

#' Create plot-aware stratified split
create_plot_aware_split <- function(data, train_prop = 0.70, val_prop = 0.15, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Get unique plots with stratification variables
  plot_info <- data %>%
    group_by(Plot_ID) %>%
    summarise(
      yield = first(yield),
      Spacing = first(Spacing),
      Leaf_type = first(Leaf_type),
      .groups = "drop"
    )
  
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
  
  unique_plots <- unique(data$Plot_ID)
  n_plots <- length(unique_plots)
  
  sampled_plots <- sample(unique_plots, size = n_plots, replace = TRUE)
  
  boot_data <- data.frame()
  for (plot in sampled_plots) {
    plot_data <- data %>% filter(Plot_ID == plot)
    boot_data <- rbind(boot_data, plot_data)
  }
  
  return(boot_data)
}

################################################################################
# SECTION 4: MODEL TRAINING FUNCTIONS
################################################################################

#' Train single model
train_single_model <- function(model_name, X_train, y_train, X_val, y_val, 
                               X_test, y_test, feature_names) {
  
  tryCatch({
    
    # Ensure no NA values in data
    if (any(is.na(X_train)) || any(is.na(y_train))) {
      stop("NA values detected in training data")
    }
    
    if (model_name == "lasso") {
      cv_model <- cv.glmnet(X_train, y_train, alpha = 1, nfolds = 5)
      model <- glmnet(X_train, y_train, alpha = 1, lambda = cv_model$lambda.min)
      pred_test <- predict(model, X_test)[, 1]
      
      coefs <- coef(model)[-1, 1]
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
      model <- randomForest(X_train, y_train, ntree = 500, 
                            mtry = sqrt(ncol(X_train)),
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
      
      imp_matrix <- xgb.importance(model = model)
      importance <- setNames(imp_matrix$Gain, imp_matrix$Feature)
      
    } else if (model_name == "svm") {
      model <- svm(X_train, y_train, kernel = "radial", 
                   gamma = 1/ncol(X_train))
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
# SECTION 5: MAIN TRAINING WITH CONFIDENCE INTERVALS
################################################################################

#' Train models with bootstrap confidence intervals
train_models_with_ci <- function(data, dataset_name, target_info, 
                                 n_reps = N_REPETITIONS, n_bootstrap = N_BOOTSTRAP) {
  
  log_message(sprintf("Training models for %s - %s with %d repetitions", 
                      dataset_name, target_info, n_reps))
  
  if (nrow(data) < 30) {
    log_message(sprintf("  ERROR: Insufficient data for modeling (%d rows)", nrow(data)), "ERROR")
    return(NULL)
  }
  
  # Define models
  models <- c("lasso", "ridge", "lm", "pls", "rf", "xgboost", "svm")
  
  # Progress tracking
  pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
  
  # Sequential execution for cluster stability
  rep_results <- list()
  
  for (rep in 1:n_reps) {
    setTxtProgressBar(pb, rep)
    
    rep_result <- tryCatch({
      # Create plot-aware split
      splits <- create_plot_aware_split(data, seed = rep * 1000)
      
      if (nrow(splits$train) < 20) {
        return(list(error = sprintf("Training set too small: %d rows", nrow(splits$train))))
      }
      
      # Prepare features
      all_cols <- names(splits$train)
      exclude_cols <- c("Plot_ID", "Spacing", "Leaf_type", "yield", "DAP", 
                        "n_timepoints", "mean_dap")
      feature_cols <- setdiff(all_cols, exclude_cols)
      feature_cols <- feature_cols[sapply(splits$train[feature_cols], is.numeric)]
      
      if (length(feature_cols) < 5) {
        return(list(error = sprintf("Too few features: %d", length(feature_cols))))
      }
      
      X_train <- as.matrix(splits$train[, feature_cols])
      y_train <- splits$train$yield
      X_val <- as.matrix(splits$val[, feature_cols])
      y_val <- splits$val$yield
      X_test <- as.matrix(splits$test[, feature_cols])
      y_test <- splits$test$yield
      
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
          
          # Fixed: Ensure all vectors have same length
          pred_df <- data.frame(
            Plot_ID = splits$test$Plot_ID,
            Observed = y_test,
            Predicted = model_results$predictions,
            Model = model_name,
            Repetition = rep,
            stringsAsFactors = FALSE
          )
          
          if (nrow(pred_df) == length(y_test)) {
            rep_predictions[[model_name]] <- pred_df
          }
          
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
      return(list(error = e$message))
    })
    
    rep_results[[rep]] <- rep_result
  }
  
  close(pb)
  
  # Check for errors
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
    return(NULL)
  }
  
  # Get feature columns
  all_cols <- names(data)
  exclude_cols <- c("Plot_ID", "Spacing", "Leaf_type", "yield", "DAP", 
                    "n_timepoints", "mean_dap")
  feature_cols <- setdiff(all_cols, exclude_cols)
  feature_cols <- feature_cols[sapply(data[feature_cols], is.numeric)]
  
  # Aggregate results
  aggregated_results <- aggregate_repetition_results(
    rep_results, models, dataset_name, target_info, n_successful
  )
  
  # Bootstrap feature importance
  if (n_bootstrap > 0 && length(feature_cols) > 0) {
    log_message("  Calculating bootstrap feature importance...")
    aggregated_results$bootstrap_importance <- calculate_bootstrap_importance(
      data, feature_cols, models = c("rf", "xgboost", "lasso", "ridge"),
      n_bootstrap = n_bootstrap
    )
  }
  
  return(aggregated_results)
}

#' Aggregate results across repetitions
aggregate_repetition_results <- function(rep_results, models, dataset_name, 
                                         target_info, n_reps) {
  
  aggregated <- list()
  
  for (model_name in models) {
    model_metrics <- list()
    
    for (rep in 1:length(rep_results)) {
      if (!is.null(rep_results[[rep]]$metrics[[model_name]])) {
        model_metrics[[rep]] <- rep_results[[rep]]$metrics[[model_name]]
      }
    }
    
    if (length(model_metrics) > 0) {
      metrics_df <- do.call(rbind, lapply(model_metrics, as.data.frame))
      
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
        target_info = target_info
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

#' Calculate bootstrap feature importance
calculate_bootstrap_importance <- function(data, feature_cols, models, n_bootstrap) {
  
  importance_results <- list()
  
  # Ensure data has no NA values
  data_clean <- data %>%
    select(Plot_ID, all_of(feature_cols), yield) %>%
    na.omit()
  
  if (nrow(data_clean) < 30) {
    log_message("    Insufficient clean data for bootstrap importance", "WARNING")
    return(NULL)
  }
  
  feature_cols_clean <- intersect(feature_cols, names(data_clean))
  
  for (model_name in models) {
    if (model_name %in% c("rf", "xgboost", "lasso", "ridge")) {
      
      # Sequential bootstrap for cluster
      boot_importance <- matrix(NA, nrow = n_bootstrap, ncol = length(feature_cols_clean))
      colnames(boot_importance) <- feature_cols_clean
      
      pb_boot <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
      
      for (b in 1:n_bootstrap) {
        setTxtProgressBar(pb_boot, b)
        
        tryCatch({
          boot_data <- bootstrap_plot_sample(data_clean)
          
          X_boot <- as.matrix(boot_data[, feature_cols_clean])
          y_boot <- boot_data$yield
          
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
            full_importance <- rep(0, length(feature_cols_clean))
            names(full_importance) <- feature_cols_clean
            full_importance[names(importance)] <- importance
            importance <- full_importance
            
          } else if (model_name %in% c("lasso", "ridge")) {
            alpha <- ifelse(model_name == "lasso", 1, 0)
            cv_model <- cv.glmnet(X_boot, y_boot, alpha = alpha, nfolds = 5)
            model <- glmnet(X_boot, y_boot, alpha = alpha, lambda = cv_model$lambda.min)
            coefs <- coef(model)[-1, 1]
            importance <- abs(coefs)
            names(importance) <- feature_cols_clean
          }
          
          boot_importance[b, ] <- importance
          
        }, error = function(e) {
          # Skip this bootstrap iteration
          warning(sprintf("Bootstrap iteration %d failed for %s: %s", b, model_name, e$message))
        })
      }
      
      close(pb_boot)
      
      if (!all(is.na(boot_importance))) {
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
        
        importance_results[[model_name]] <- imp_summary
      }
    }
  }
  
  return(importance_results)
}

################################################################################
# SECTION 6: STATISTICAL TESTING
################################################################################

#' Perform pairwise statistical tests between models
compare_models_statistically <- function(results_list, metric = "R2") {
  
  model_names <- names(results_list)
  n_models <- length(model_names)
  
  p_matrix <- matrix(1, nrow = n_models, ncol = n_models)
  rownames(p_matrix) <- colnames(p_matrix) <- model_names
  
  metric_values <- list()
  for (model in model_names) {
    if (!is.null(results_list[[model]]$all_predictions)) {
      predictions <- results_list[[model]]$all_predictions
      
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
  
  for (i in 1:(n_models-1)) {
    for (j in (i+1):n_models) {
      if (length(metric_values[[i]]) > 0 && length(metric_values[[j]]) > 0) {
        test_result <- wilcox.test(metric_values[[i]], metric_values[[j]], 
                                   paired = FALSE, exact = FALSE)
        p_matrix[i, j] <- p_matrix[j, i] <- test_result$p.value
      }
    }
  }
  
  return(p_matrix)
}

################################################################################
# SECTION 7: MAIN ANALYSIS - INDIVIDUAL TIMEPOINTS
################################################################################

log_message("\n=== YIELD ANALYSIS WITH ENHANCED STATISTICAL RIGOR ===")

# Process individual timepoints
all_timepoint_results <- list()
all_timepoint_bootstrap <- list()

# Create checkpoint file to track progress
timepoint_checkpoint_file <- file.path(output_dir, "checkpoints", "timepoint_checkpoint.rds")
if (file.exists(timepoint_checkpoint_file)) {
  log_message("Loading timepoint checkpoint data...")
  checkpoint_data <- readRDS(timepoint_checkpoint_file)
  all_timepoint_results <- checkpoint_data$results
  all_timepoint_bootstrap <- checkpoint_data$bootstrap
  completed_analyses <- checkpoint_data$completed
} else {
  completed_analyses <- character()
}

for (dataset_name in names(loaded_data)) {
  log_message(sprintf("\nProcessing timepoints for %s", dataset_name))
  
  data <- loaded_data[[dataset_name]]
  unique_timepoints <- unique(data$Timepoint)
  
  log_message(sprintf("  Found %d timepoints", length(unique_timepoints)))
  
  dataset_timepoint_results <- list()
  
  for (tp in unique_timepoints) {
    analysis_key <- paste(dataset_name, tp, sep = "_")
    
    # Skip if already completed
    if (analysis_key %in% completed_analyses) {
      log_message(sprintf("\n  Skipping timepoint %s (already completed)", tp))
      if (!is.null(all_timepoint_results[[dataset_name]][[tp]])) {
        dataset_timepoint_results[[tp]] <- all_timepoint_results[[dataset_name]][[tp]]
      }
      next
    }
    
    log_message(sprintf("\n  Processing timepoint %s", tp))
    
    # Prepare data
    tp_data <- tryCatch({
      prepare_yield_timepoint_data(data, tp)
    }, error = function(e) {
      log_message(sprintf("    ERROR preparing data: %s", e$message), "ERROR")
      NULL
    })
    
    if (!is.null(tp_data) && nrow(tp_data) > 50) {
      # Get mean DAP
      mean_dap <- mean(tp_data$DAP, na.rm = TRUE)
      
      # Train models with CI
      model_results <- tryCatch({
        train_models_with_ci(
          tp_data, 
          dataset_name,
          target_info = tp,
          n_reps = N_REPETITIONS,
          n_bootstrap = N_BOOTSTRAP
        )
      }, error = function(e) {
        log_message(sprintf("    ERROR in modeling: %s", e$message), "ERROR")
        NULL
      })
      
      if (!is.null(model_results)) {
        model_results$mean_dap <- mean_dap
        dataset_timepoint_results[[tp]] <- model_results
        
        if (!is.null(model_results$bootstrap_importance)) {
          all_timepoint_bootstrap[[paste(dataset_name, tp, sep = "_")]] <- 
            model_results$bootstrap_importance
        }
        
        # Mark as completed and save checkpoint
        completed_analyses <- c(completed_analyses, analysis_key)
        
        # Update results
        all_timepoint_results[[dataset_name]] <- dataset_timepoint_results
        
        # Save checkpoint
        saveRDS(list(
          results = all_timepoint_results,
          bootstrap = all_timepoint_bootstrap,
          completed = completed_analyses
        ), timepoint_checkpoint_file)
        
        # Also save individual result file
        saveRDS(model_results, 
                file.path(output_dir, "individual_timepoints", 
                          sprintf("%s_%s_results.rds", dataset_name, tp)))
        
        log_message(sprintf("    Saved results for %s", tp))
      }
    } else {
      log_message(sprintf("    Insufficient data for %s", tp), "WARNING")
    }
  }
  
  all_timepoint_results[[dataset_name]] <- dataset_timepoint_results
  
  # Save intermediate results after each dataset
  saveRDS(all_timepoint_results, file.path(output_dir, "models", 
                                            "yield_timepoint_intermediate_results.rds"))
  log_message(sprintf("  Saved intermediate results for %s", dataset_name))
}

# Save final timepoint results
saveRDS(all_timepoint_results, file.path(output_dir, "models", 
                                          "yield_timepoint_results_with_ci.rds"))

################################################################################
# SECTION 8: MAIN ANALYSIS - TEMPORAL AGGREGATION
################################################################################

log_message("\n=== TEMPORAL AGGREGATION ANALYSIS ===")

all_temporal_results <- list()
all_temporal_bootstrap <- list()

# Create checkpoint file for temporal analysis
temporal_checkpoint_file <- file.path(output_dir, "checkpoints", "temporal_checkpoint.rds")
if (file.exists(temporal_checkpoint_file)) {
  log_message("Loading temporal checkpoint data...")
  temporal_checkpoint <- readRDS(temporal_checkpoint_file)
  all_temporal_results <- temporal_checkpoint$results
  all_temporal_bootstrap <- temporal_checkpoint$bootstrap
  completed_temporal_analyses <- temporal_checkpoint$completed
} else {
  completed_temporal_analyses <- character()
}

for (dataset_name in names(loaded_data)) {
  # Skip if already completed
  if (dataset_name %in% completed_temporal_analyses) {
    log_message(sprintf("\nSkipping temporal aggregation for %s (already completed)", dataset_name))
    next
  }
  
  log_message(sprintf("\nProcessing temporal aggregation for %s", dataset_name))
  
  # Prepare temporal aggregation data
  temporal_data <- tryCatch({
    prepare_yield_temporal_aggregation(loaded_data[[dataset_name]])
  }, error = function(e) {
    log_message(sprintf("  ERROR preparing temporal data: %s", e$message), "ERROR")
    NULL
  })
  
  if (!is.null(temporal_data) && nrow(temporal_data) > 50) {
    # Train models with CI
    model_results <- tryCatch({
      train_models_with_ci(
        temporal_data, 
        dataset_name,
        target_info = "temporal_aggregation",
        n_reps = N_REPETITIONS,
        n_bootstrap = N_BOOTSTRAP
      )
    }, error = function(e) {
      log_message(sprintf("  ERROR in temporal modeling: %s", e$message), "ERROR")
      NULL
    })
    
    if (!is.null(model_results)) {
      all_temporal_results[[dataset_name]] <- model_results
      
      if (!is.null(model_results$bootstrap_importance)) {
        all_temporal_bootstrap[[dataset_name]] <- model_results$bootstrap_importance
      }
      
      # Mark as completed and save checkpoint
      completed_temporal_analyses <- c(completed_temporal_analyses, dataset_name)
      
      # Save checkpoint
      saveRDS(list(
        results = all_temporal_results,
        bootstrap = all_temporal_bootstrap,
        completed = completed_temporal_analyses
      ), temporal_checkpoint_file)
      
      # Save individual temporal result
      saveRDS(model_results, 
              file.path(output_dir, "temporal_analysis", 
                        sprintf("%s_temporal_results.rds", dataset_name)))
      
      log_message(sprintf("  Saved temporal results for %s", dataset_name))
    }
  } else {
    log_message(sprintf("  Insufficient data for temporal aggregation in %s", dataset_name), "WARNING")
  }
  
  # Save intermediate results after each dataset
  saveRDS(all_temporal_results, file.path(output_dir, "models", 
                                           "yield_temporal_intermediate_results.rds"))
}

# Save final temporal results
saveRDS(all_temporal_results, file.path(output_dir, "models", 
                                         "yield_temporal_results_with_ci.rds"))

# Combine and save all bootstrap results
all_bootstrap_combined <- c(all_timepoint_bootstrap, all_temporal_bootstrap)
if (length(all_bootstrap_combined) > 0) {
  saveRDS(all_bootstrap_combined, 
          file.path(output_dir, "bootstrap_results", "feature_importance_bootstrap.rds"))
}

################################################################################
# SECTION 9: COMPILE PERFORMANCE SUMMARIES
################################################################################

log_message("\n=== Compiling Performance Summaries ===")

# Wrap the entire summary section in tryCatch
tryCatch({
  
  # Compile timepoint results
  timepoint_summary_ci <- data.frame()
  
  for (dataset_name in names(all_timepoint_results)) {
    dataset_results <- all_timepoint_results[[dataset_name]]
    
    for (tp in names(dataset_results)) {
      tp_result <- dataset_results[[tp]]
      mean_dap <- tp_result$mean_dap
      
      for (model_name in names(tp_result)) {
        if (model_name %in% c("all_predictions", "n_repetitions", "bootstrap_importance", "mean_dap")) {
          next
        }
        
        model_result <- tp_result[[model_name]]
        
        if (!is.null(model_result$metrics)) {
          metrics <- model_result$metrics
          
          timepoint_summary_ci <- rbind(timepoint_summary_ci, data.frame(
            Dataset = dataset_name,
            Timepoint = tp,
            Mean_DAP = mean_dap,
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
  if (nrow(timepoint_summary_ci) > 0) {
    timepoint_summary_ci <- timepoint_summary_ci %>%
      mutate(
        R2_formatted = sprintf("%.3f [%.3f, %.3f]", R2_mean, R2_ci_lower, R2_ci_upper),
        RMSE_formatted = sprintf("%.3f [%.3f, %.3f]", RMSE_mean, RMSE_ci_lower, RMSE_ci_upper),
        MAE_formatted = sprintf("%.3f [%.3f, %.3f]", MAE_mean, MAE_ci_lower, MAE_ci_upper)
      )
    
    # Save timepoint summary
    write_csv(timepoint_summary_ci, 
              file.path(output_dir, "confidence_intervals", 
                        "yield_timepoint_performance_with_ci.csv"))
    log_message("  Saved timepoint performance summary")
  }
  
  # Compile temporal aggregation results
  temporal_summary_ci <- data.frame()
  
  for (dataset_name in names(all_temporal_results)) {
    dataset_result <- all_temporal_results[[dataset_name]]
    
    for (model_name in names(dataset_result)) {
      if (model_name %in% c("all_predictions", "n_repetitions", "bootstrap_importance")) {
        next
      }
      
      model_result <- dataset_result[[model_name]]
      
      if (!is.null(model_result$metrics)) {
        metrics <- model_result$metrics
        
        temporal_summary_ci <- rbind(temporal_summary_ci, data.frame(
          Dataset = dataset_name,
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
  
  # Add formatted columns
  if (nrow(temporal_summary_ci) > 0) {
    temporal_summary_ci <- temporal_summary_ci %>%
      mutate(
        R2_formatted = sprintf("%.3f [%.3f, %.3f]", R2_mean, R2_ci_lower, R2_ci_upper),
        RMSE_formatted = sprintf("%.3f [%.3f, %.3f]", RMSE_mean, RMSE_ci_lower, RMSE_ci_upper),
        MAE_formatted = sprintf("%.3f [%.3f, %.3f]", MAE_mean, MAE_ci_lower, MAE_ci_upper)
      )
    
    # Save temporal summary
    write_csv(temporal_summary_ci, 
              file.path(output_dir, "confidence_intervals", 
                        "yield_temporal_performance_with_ci.csv"))
    log_message("  Saved temporal performance summary")
  }
  
}, error = function(e) {
  log_message(sprintf("ERROR in performance summary compilation: %s", e$message), "ERROR")
  log_message("Continuing with analysis despite summary error...", "WARNING")
})

################################################################################
# SECTION 10: STATISTICAL COMPARISONS
################################################################################

log_message("\n=== Performing Statistical Comparisons ===")

# Wrap in tryCatch to prevent crashes
tryCatch({
  
  # Timepoint comparisons
  timepoint_comparisons <- list()
  
  for (dataset_name in names(all_timepoint_results)) {
    dataset_results <- all_timepoint_results[[dataset_name]]
    
    for (tp in names(dataset_results)) {
      tp_result <- dataset_results[[tp]]
      
      model_names <- setdiff(names(tp_result), 
                             c("all_predictions", "n_repetitions", "bootstrap_importance", "mean_dap"))
      
      if (length(model_names) > 1) {
        model_list <- list()
        for (model in model_names) {
          model_list[[model]] <- tp_result[[model]]
        }
        
        tryCatch({
          p_matrix <- compare_models_statistically(model_list, metric = "R2")
          
          comparison_key <- paste(dataset_name, tp, sep = "_")
          timepoint_comparisons[[comparison_key]] <- p_matrix
          
          write.csv(p_matrix, 
                    file.path(output_dir, "statistical_tests",
                              sprintf("pvalues_%s_%s.csv", dataset_name, tp)))
        }, error = function(e) {
          log_message(sprintf("  Error in statistical comparison for %s %s: %s", 
                              dataset_name, tp, e$message), "WARNING")
        })
      }
    }
  }
  
  # Temporal aggregation comparisons
  temporal_comparisons <- list()
  
  for (dataset_name in names(all_temporal_results)) {
    dataset_result <- all_temporal_results[[dataset_name]]
    
    model_names <- setdiff(names(dataset_result), 
                           c("all_predictions", "n_repetitions", "bootstrap_importance"))
    
    if (length(model_names) > 1) {
      model_list <- list()
      for (model in model_names) {
        model_list[[model]] <- dataset_result[[model]]
      }
      
      tryCatch({
        p_matrix <- compare_models_statistically(model_list, metric = "R2")
        
        temporal_comparisons[[dataset_name]] <- p_matrix
        
        write.csv(p_matrix, 
                  file.path(output_dir, "statistical_tests",
                            sprintf("pvalues_%s_temporal.csv", dataset_name)))
      }, error = function(e) {
        log_message(sprintf("  Error in statistical comparison for %s temporal: %s", 
                            dataset_name, e$message), "WARNING")
      })
    }
  }
  
  log_message("  Statistical comparisons completed (where possible)")
  
}, error = function(e) {
  log_message(sprintf("ERROR in statistical comparisons: %s", e$message), "ERROR")
  log_message("Continuing with analysis despite comparison error...", "WARNING")
})

################################################################################
# SECTION 11: VISUALIZATIONS
################################################################################

log_message("\n=== Creating Visualizations ===")

# Wrap visualization section in tryCatch
tryCatch({
  
  # Performance comparison with error bars
  create_performance_comparison_plot <- function(performance_data, title_suffix = "") {
    
    p <- ggplot(performance_data, aes(x = Model, y = R2_mean, color = Dataset)) +
      geom_point(position = position_dodge(0.5), size = 3) +
      geom_errorbar(aes(ymin = R2_ci_lower, ymax = R2_ci_upper), 
                    width = 0.2, position = position_dodge(0.5)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      labs(title = paste0("Model Performance with 95% Confidence Intervals", title_suffix),
           x = "Model", 
           y = "R² (mean with 95% CI)") +
      scale_color_brewer(palette = "Set2")
    
    return(p)
  }
  
  # Feature importance with uncertainty
  create_feature_importance_plot <- function(importance_data, model_name, top_n = 20) {
    
    top_features <- importance_data %>%
      arrange(desc(importance_mean)) %>%
      slice_head(n = top_n) %>%
      mutate(feature = factor(feature, levels = rev(feature)))
    
    p <- ggplot(top_features, aes(x = importance_mean, y = feature)) +
      geom_point(size = 3) +
      geom_errorbarh(aes(xmin = importance_ci_lower, xmax = importance_ci_upper), 
                     height = 0.2) +
      theme_minimal() +
      labs(title = sprintf("Feature Importance with 95%% CI - %s", model_name),
           x = "Importance (mean with 95% CI)", 
           y = "Feature") +
      theme(axis.text.y = element_text(size = 8))
    
    return(p)
  }
  
  # Temporal trend plot
  create_temporal_trend_plot <- function(timepoint_data, dataset_name) {
    
    p <- ggplot(timepoint_data %>% filter(Dataset == dataset_name), 
                aes(x = Mean_DAP, y = R2_mean, color = Model)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      geom_ribbon(aes(ymin = R2_ci_lower, ymax = R2_ci_upper, fill = Model), 
                  alpha = 0.2) +
      theme_minimal() +
      labs(title = sprintf("Yield Prediction Performance Over Time - %s", dataset_name),
           x = "Days After Planting",
           y = "R² (mean with 95% CI)") +
      scale_color_brewer(palette = "Set2") +
      scale_fill_brewer(palette = "Set2")
    
    return(p)
  }
  
  # Create visualizations
  pdf(file.path(output_dir, "visualizations", "yield_analysis_enhanced.pdf"), 
      width = 14, height = 10)
  
  # Title page
  plot.new()
  text(0.5, 0.7, "Yield Prediction Analysis with Statistical Rigor", cex = 2, font = 2)
  text(0.5, 0.5, sprintf("%d Bootstrap Iterations, %d Repeated Holdouts", 
                         N_BOOTSTRAP, N_REPETITIONS), cex = 1.5)
  text(0.5, 0.3, format(Sys.time(), "%Y-%m-%d"), cex = 1.2)
  
  # Temporal aggregation performance
  if (exists("temporal_summary_ci") && nrow(temporal_summary_ci) > 0) {
    tryCatch({
      p1 <- create_performance_comparison_plot(temporal_summary_ci, 
                                               " - Temporal Aggregation")
      print(p1)
    }, error = function(e) {
      log_message(sprintf("  Error creating temporal plot: %s", e$message), "WARNING")
    })
  }
  
  # Timepoint performance by dataset
  if (exists("timepoint_summary_ci") && nrow(timepoint_summary_ci) > 0) {
    for (dataset in unique(timepoint_summary_ci$Dataset)) {
      tryCatch({
        p_trend <- create_temporal_trend_plot(timepoint_summary_ci, dataset)
        print(p_trend)
      }, error = function(e) {
        log_message(sprintf("  Error creating trend plot for %s: %s", dataset, e$message), "WARNING")
      })
    }
  }
  
  # Feature importance for each model with bootstrap CI
  for (dataset in names(all_temporal_bootstrap)) {
    for (model in names(all_temporal_bootstrap[[dataset]])) {
      tryCatch({
        imp_data <- all_temporal_bootstrap[[dataset]][[model]]
        
        if (!is.null(imp_data) && nrow(imp_data) > 0) {
          p <- create_feature_importance_plot(imp_data, 
                                              sprintf("%s - %s", dataset, model))
          print(p)
        }
      }, error = function(e) {
        log_message(sprintf("  Error creating importance plot for %s %s: %s", 
                            dataset, model, e$message), "WARNING")
      })
    }
  }
  
  # Best timepoint analysis
  if (exists("timepoint_summary_ci") && nrow(timepoint_summary_ci) > 0) {
    tryCatch({
      best_timepoints <- timepoint_summary_ci %>%
        group_by(Dataset, Model) %>%
        slice_max(R2_mean, n = 1) %>%
        ungroup()
      
      p_best <- ggplot(best_timepoints, aes(x = Model, y = Mean_DAP, fill = Dataset)) +
        geom_bar(stat = "identity", position = "dodge") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Optimal Timepoint (DAP) for Yield Prediction by Model",
             x = "Model", y = "Days After Planting") +
        scale_fill_brewer(palette = "Set1")
      
      print(p_best)
    }, error = function(e) {
      log_message(sprintf("  Error creating best timepoint plot: %s", e$message), "WARNING")
    })
  }
  
  # Comparison: Timepoint vs Temporal
  if (exists("timepoint_summary_ci") && nrow(timepoint_summary_ci) > 0 && 
      exists("temporal_summary_ci") && nrow(temporal_summary_ci) > 0) {
    tryCatch({
      best_timepoint <- timepoint_summary_ci %>%
        group_by(Dataset, Model) %>%
        slice_max(R2_mean, n = 1) %>%
        select(Dataset, Model, R2_mean, R2_ci_lower, R2_ci_upper) %>%
        rename(Timepoint_R2 = R2_mean, 
               Timepoint_Lower = R2_ci_lower,
               Timepoint_Upper = R2_ci_upper)
      
      comparison <- temporal_summary_ci %>%
        select(Dataset, Model, R2_mean, R2_ci_lower, R2_ci_upper) %>%
        rename(Temporal_R2 = R2_mean,
               Temporal_Lower = R2_ci_lower,
               Temporal_Upper = R2_ci_upper) %>%
        inner_join(best_timepoint, by = c("Dataset", "Model"))
      
      p_comp <- ggplot(comparison, aes(x = Timepoint_R2, y = Temporal_R2, 
                                       color = Dataset, shape = Model)) +
        geom_point(size = 3) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        geom_errorbar(aes(ymin = Temporal_Lower, ymax = Temporal_Upper), 
                      width = 0.02, alpha = 0.5) +
        geom_errorbarh(aes(xmin = Timepoint_Lower, xmax = Timepoint_Upper), 
                       height = 0.02, alpha = 0.5) +
        theme_minimal() +
        labs(title = "Yield Prediction: Temporal Aggregation vs Best Individual Timepoint",
             subtitle = "With 95% confidence intervals",
             x = "Best Timepoint R²", y = "Temporal Aggregation R²") +
        coord_fixed() +
        xlim(0, 1) + ylim(0, 1)
      
      print(p_comp)
    }, error = function(e) {
      log_message(sprintf("  Error creating comparison plot: %s", e$message), "WARNING")
    })
  }
  
  dev.off()
  log_message("  Visualizations completed")
  
}, error = function(e) {
  log_message(sprintf("ERROR in visualization creation: %s", e$message), "ERROR")
  log_message("Continuing with analysis despite visualization error...", "WARNING")
  try(dev.off(), silent = TRUE)  # Close any open device
})

################################################################################
# SECTION 12: SUMMARY REPORT
################################################################################

log_message("\n=== Creating Summary Report ===")

# Wrap summary report in tryCatch
tryCatch({
  
  sink(file.path(output_dir, "yield_analysis_enhanced_summary.txt"))
  
  cat("YIELD PREDICTION ANALYSIS WITH ENHANCED STATISTICAL RIGOR\n")
  cat("=========================================================\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("ANALYSIS CONFIGURATION:\n")
  cat(sprintf("- Repeated holdout validation: %d repetitions\n", N_REPETITIONS))
  cat(sprintf("- Bootstrap iterations for feature importance: %d\n", N_BOOTSTRAP))
  cat(sprintf("- Confidence level: %.0f%%\n", CONFIDENCE_LEVEL * 100))
  cat("- Train/Val/Test split: 70%/15%/15%\n")
  cat("- Plot-aware sampling to prevent data leakage\n\n")
  
  cat("KEY FINDINGS:\n\n")
  
  # Best models for temporal aggregation
  if (exists("temporal_summary_ci") && nrow(temporal_summary_ci) > 0) {
    cat("1. BEST MODELS FOR TEMPORAL AGGREGATION (with 95% CI):\n")
    best_temporal <- temporal_summary_ci %>%
      group_by(Dataset) %>%
      slice_max(R2_mean, n = 1) %>%
      select(Dataset, Model, R2_formatted, RMSE_formatted)
    
    for (i in 1:nrow(best_temporal)) {
      cat(sprintf("   %s: %s model, R² = %s, RMSE = %s\n",
                  best_temporal$Dataset[i], 
                  best_temporal$Model[i],
                  best_temporal$R2_formatted[i],
                  best_temporal$RMSE_formatted[i]))
    }
  }
  
  # Best timepoints
  if (exists("timepoint_summary_ci") && nrow(timepoint_summary_ci) > 0) {
    cat("\n2. OPTIMAL TIMEPOINTS FOR YIELD PREDICTION:\n")
    optimal_timepoints <- timepoint_summary_ci %>%
      group_by(Dataset) %>%
      slice_max(R2_mean, n = 1) %>%
      select(Dataset, Timepoint, Mean_DAP, Model, R2_formatted)
    
    for (i in 1:nrow(optimal_timepoints)) {
      cat(sprintf("   %s: %s (DAP %.0f), %s model, R² = %s\n",
                  optimal_timepoints$Dataset[i],
                  optimal_timepoints$Timepoint[i],
                  optimal_timepoints$Mean_DAP[i],
                  optimal_timepoints$Model[i],
                  optimal_timepoints$R2_formatted[i]))
    }
  }
  
  cat("\n3. STATISTICAL SIGNIFICANCE:\n")
  # Count significant differences
  n_sig_timepoint <- 0
  if (exists("timepoint_comparisons")) {
    for (comparison in names(timepoint_comparisons)) {
      p_matrix <- timepoint_comparisons[[comparison]]
      n_sig <- sum(p_matrix < 0.05 & p_matrix != 1, na.rm = TRUE) / 2
      n_sig_timepoint <- n_sig_timepoint + n_sig
    }
  }
  
  n_sig_temporal <- 0
  if (exists("temporal_comparisons")) {
    for (comparison in names(temporal_comparisons)) {
      p_matrix <- temporal_comparisons[[comparison]]
      n_sig <- sum(p_matrix < 0.05 & p_matrix != 1, na.rm = TRUE) / 2
      n_sig_temporal <- n_sig_temporal + n_sig
    }
  }
  
  cat(sprintf("   Timepoint analyses: %d significant model differences\n", n_sig_timepoint))
  cat(sprintf("   Temporal aggregation: %d significant model differences\n", n_sig_temporal))
  
  cat("\n4. MODEL CONSISTENCY:\n")
  # All results combined
  all_results <- data.frame()
  if (exists("timepoint_summary_ci") && nrow(timepoint_summary_ci) > 0) {
    all_results <- rbind(all_results, timepoint_summary_ci %>% mutate(Type = "Timepoint"))
  }
  if (exists("temporal_summary_ci") && nrow(temporal_summary_ci) > 0) {
    all_results <- rbind(all_results, temporal_summary_ci %>% mutate(Type = "Temporal"))
  }
  
  if (nrow(all_results) > 0) {
    model_consistency <- all_results %>%
      group_by(Model) %>%
      summarise(
        Mean_R2 = mean(R2_mean, na.rm = TRUE),
        CV_R2 = sd(R2_mean, na.rm = TRUE) / mean(R2_mean, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(CV_R2)
    
    cat("   Models ranked by consistency (lower CV = more consistent):\n")
    for (i in 1:nrow(model_consistency)) {
      cat(sprintf("   %d. %s: Mean R² = %.3f, CV = %.3f\n",
                  i,
                  model_consistency$Model[i],
                  model_consistency$Mean_R2[i],
                  model_consistency$CV_R2[i]))
    }
  }
  
  cat("\n5. DATA PROCESSING SUMMARY:\n")
  
  # Count successful analyses
  n_timepoint_success <- 0
  n_timepoint_total <- 0
  for (dataset in names(all_timepoint_results)) {
    n_timepoint_total <- n_timepoint_total + length(unique(loaded_data[[dataset]]$Timepoint))
    n_timepoint_success <- n_timepoint_success + length(all_timepoint_results[[dataset]])
  }
  
  n_temporal_success <- length(all_temporal_results)
  n_temporal_total <- length(loaded_data)
  
  cat(sprintf("   Individual timepoints: %d/%d successfully analyzed\n", 
              n_timepoint_success, n_timepoint_total))
  cat(sprintf("   Temporal aggregation: %d/%d datasets successfully analyzed\n", 
              n_temporal_success, n_temporal_total))
  
  sink()
  
  log_message("  Summary report created successfully")
  
}, error = function(e) {
  log_message(sprintf("ERROR creating summary report: %s", e$message), "ERROR")
  try(sink(), silent = TRUE)  # Close sink if open
})

# Clean up (parallel cluster already disabled)
log_message("\n=== Phase 2 Yield Analysis Complete (Enhanced Version) ===")
log_message(sprintf("Results saved to: %s", output_dir))
log_message("All analyses include confidence intervals and statistical significance tests")
log_message("Analysis includes checkpoint recovery and incremental saving")