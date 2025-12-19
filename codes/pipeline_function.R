# ==============================================================================
# Pipeline Functions for Gene Expression Analysis
# ==============================================================================
#
# This script contains a collection of functions to perform a machine learning
# pipeline for gene expression data. The steps include data resampling, 
# Gene Set Variation Analysis (GSVA), Random Forest model training and 
# evaluation, and calibration curve generation.
#
# ==============================================================================
# I. Data Resampling
# ==============================================================================

#' Resample Gene Expression Data
#'
#' This function takes a gene expression matrix and an annotation file, and
#' performs over-sampling to balance the classes.
#'
#' @param gene_expression_file A matrix or data frame with genes in rows and samples in columns.
#' @param annot_file A data frame with sample annotations, including 'ID' and 'class' columns.
#' @param seed An optional integer for setting the random seed for reproducibility.
#'
#' @return A list containing the over-sampled expression matrix (`expr_over`) and 
#'         the corresponding annotation data frame (`annot_over`).
#'
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr left_join select mutate
#' @importFrom ROSE ovun.sample
#'
get_resample_data <- function(gene_expression_file, annot_file, seed = NULL) {
  
  # Initialize an empty list to store the results
  list_file <- list()
  
  # Combine gene expression and annotation data
  gene_annot_matrix <- gene_expression_file %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    left_join(annot_file %>% dplyr::select(ID, class))
  
  # Ensure column names are unique and syntactically valid
  colnames(gene_annot_matrix) <- make.names(colnames(gene_annot_matrix), unique = TRUE)
  gene_annot_matrix$class <- as.factor(gene_annot_matrix$class)
  
  # Define a helper function for resampling
  resample_data <- function(data, method) {
    ovun.sample(class ~ ., data = data, method = method)$data
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Perform over-sampling
  data_over <- resample_data(data = gene_annot_matrix, method = "over")
  
  # Prepare the over-sampled expression matrix
  list_file$expr_over <- data_over %>%
    mutate(new_id = paste0("id_", 1:min(dim(data_over)))) %>%
    column_to_rownames("new_id") %>%
    dplyr::select(-class, -ID) %>%
    as.data.frame() %>%
    t()
  
  # Prepare the corresponding annotation file
  list_file$annot_over <- data.frame(ID = colnames(list_file$expr_over), class = data_over$class) %>%
    mutate(class = factor(class))
  
  return(list_file)
}


# ==============================================================================
# II. Gene Set Variation Analysis (GSVA)
# ==============================================================================

#' Perform Gene Set Variation Analysis (GSVA)
#'
#' This function calculates GSVA scores for a given gene expression matrix and a 
#' list of gene sets.
#'
#' @param tpm_matrix A matrix of TPM (Transcripts Per Million) values with genes in rows and samples in columns.
#' @param gsva_feature_list A list of character vectors, where each vector is a gene set.
#' @param annot_matrix A data frame with sample annotations.
#' @param disease_name A character string specifying the name of the disease class.
#' @param control_name A character string specifying the name of the control class.
#'
#' @return A list containing the raw GSVA scores (`gsva_score`) and a data frame 
#'         of feature scores combined with sample annotations (`feature_score`).
#'
#' @importFrom GSVA gsva
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join relocate mutate
#'
do_gsva <- function(tpm_matrix, gsva_feature_list, annot_matrix, disease_name, control_name) {
  
  list_name <- list()
  
  # Calculate GSVA scores
  list_name$gsva_score <- gsva(
    expr = tpm_matrix %>% as.matrix(),
    gset.idx.list = gsva_feature_list,
    method = "gsva",
    kcdf = "Gaussian", # Using Gaussian kernel for continuous data
    parallel.sz = 8L
  )
  
  # Combine GSVA scores with sample annotations
  list_name$feature_score <- list_name$gsva_score %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>%
    left_join(annot_matrix %>% dplyr::select(ID, class)) %>%
    relocate(class, .after = ID) %>%
    mutate(class = factor(class, levels = c(disease_name, control_name)))
    
  return(list_name)
}


# ==============================================================================
# III. Random Forest Machine Learning Model
# ==============================================================================

#' Train and Evaluate a Random Forest Model
#'
#' This function trains a Random Forest model using `tidymodels`, performs
#' hyperparameter tuning, and evaluates the model on testing data.
#'
#' @param training_data A data frame for training the model.
#' @param testing_data_1 A data frame for the first set of testing.
#' @param testing_data_2 An optional data frame for a second set of testing.
#' @param preprocess_step A character string specifying the preprocessing method: 
#'        "PCA" (Normalization + PCA), "NOR" (Normalization only), or "NON" (None).
#' @param disease_name The name of the positive class (disease).
#' @param seed An optional integer for setting the random seed.
#'
#' @return A list containing the trained model, workflow, tuning results, 
#'         predictions, and performance metrics.
#'
#' @importFrom tidymodels recipe step_normalize step_pca step_nzv step_zv 
#'             vfold_cv rand_forest set_engine set_mode workflow add_recipe 
#'             add_model grid_max_entropy tune_grid select_best finalize_workflow
#'             fit accuracy bal_accuracy roc_auc
#' @importFrom doParallel makeCluster registerDoParallel
#' @importFrom tidyverse filter arrange select pull
#' @importFrom dplyr bind_cols
#'
do_ml_rf <- function(training_data, testing_data_1, testing_data_2, preprocess_step, disease_name, seed = NULL) {
  require(tidyverse)
  require(tidymodels)
  require(doParallel)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Set up parallel processing
  cl <- makeCluster(8)
  registerDoParallel(cl)
  
  list_name <- list()
  
  # --- 1. Preprocessing Recipe ---
  if (preprocess_step == "PCA") {
    list_name$recipe <-
      recipe(class ~ ., data = training_data[, -1]) %>%
      step_normalize(all_numeric()) %>%
      step_pca(all_numeric(), num_comp = 10) %>%
      step_nzv(all_predictors()) %>%
      step_zv(all_predictors())
  } else if (preprocess_step == "NOR") {
    list_name$recipe <-
      recipe(class ~ ., data = training_data[, -1]) %>%
      step_normalize(all_numeric()) %>%
      step_nzv(all_predictors()) %>%
      step_zv(all_predictors())
  } else if (preprocess_step == "NON") {
    list_name$recipe <-
      recipe(class ~ ., data = training_data[, -1]) %>%
      step_nzv(all_predictors()) %>%
      step_zv(all_predictors())
  }
  
  # --- 2. Cross-Validation Setup ---
  list_name$cv_folds <-
    vfold_cv(training_data[, -1], v = 10, strata = class, repeats = 3)
  
  # --- 3. Model Specification ---
  list_name$rf_model <- rand_forest(
    trees = 1000,
    mtry = tune(),
    min_n = tune()
  ) %>%
    set_engine("ranger", importance = "impurity") %>%
    set_mode("classification")
  
  # --- 4. Workflow Setup ---
  list_name$rf_workflow <-
    workflow() %>%
    add_recipe(list_name$recipe) %>%
    add_model(list_name$rf_model)
  
  # --- 5. Hyperparameter Tuning ---
  list_name$rf_grid <-
    grid_max_entropy(extract_parameter_set_dials(list_name$rf_model) %>%
                       finalize(x = training_data[, -1] %>% dplyr::select(-class)), size = 20)
  
  list_name$rf_tuned <-
    tune_grid(
      object = list_name$rf_workflow,
      resamples = list_name$cv_folds,
      grid = list_name$rf_grid,
      metrics = metric_set(accuracy, bal_accuracy),
      control = tune::control_grid(verbose = TRUE)
    )
  
  # --- 6. Finalize and Fit Model ---
  list_name$rf_best_hyper <- list_name$rf_tuned %>%
    select_best(metric = "accuracy")
  
  list_name$rf_final_workflow <-
    list_name$rf_workflow %>%
    finalize_workflow(list_name$rf_best_hyper)
  
  list_name$rf_fit <-
    list_name$rf_final_workflow %>%
    fit(data = training_data[, -1])
  
  # --- 7. Evaluation ---
  # Initialize results matrix
  list_name$result <- matrix(
    0, nrow = 1, ncol = 8,
    dimnames = list(c("rf"),
                    c("cv_accuracy", "cv_bal_accuracy",
                      "pre_accuracy_1", "pre_bal_accuracy_1", "auc_1",
                      "pre_accuracy_2", "pre_bal_accuracy_2", "auc_2"))
  )
  
  # Cross-validation results
  list_name$result[1, 1] <- (list_name$rf_tuned %>% collect_metrics() %>% filter(.metric == "accuracy") %>% arrange(-mean) %>% dplyr::select(mean))[1, ] %>% pull
  list_name$result[1, 2] <- (list_name$rf_tuned %>% collect_metrics() %>% filter(.metric == "bal_accuracy") %>% arrange(-mean) %>% dplyr::select(mean))[1, ] %>% pull
  
  # Evaluation on testing_data_1
  list_name$rf_pred_1 <-
    predict(list_name$rf_fit, testing_data_1[, -1]) %>%
    bind_cols(predict(list_name$rf_fit, testing_data_1[, -1], type = "prob")) %>%
    bind_cols(testing_data_1 %>% dplyr::select(class))
  
  list_name$result[1, 3] <- list_name$rf_pred_1 %>% accuracy(truth = class, .pred_class) %>% dplyr::select(.estimate) %>% pull
  list_name$result[1, 4] <- list_name$rf_pred_1 %>% bal_accuracy(truth = class, .pred_class) %>% dplyr::select(.estimate) %>% pull
  list_name$result[1, 5] <- list_name$rf_pred_1 %>% roc_auc(truth = class, paste0(".pred_", disease_name)) %>% dplyr::select(.estimate) %>% pull
  
  # Evaluation on testing_data_2 (if provided)
  if (!is.null(testing_data_2)) {
    list_name$rf_pred_2 <-
      predict(list_name$rf_fit, testing_data_2[, -1]) %>%
      bind_cols(predict(list_name$rf_fit, testing_data_2[, -1], type = "prob")) %>%
      bind_cols(testing_data_2 %>% dplyr::select(class))
    
    list_name$result[1, 6] <- list_name$rf_pred_2 %>% accuracy(truth = class, .pred_class) %>% dplyr::select(.estimate) %>% pull
    list_name$result[1, 7] <- list_name$rf_pred_2 %>% bal_accuracy(truth = class, .pred_class) %>% dplyr::select(.estimate) %>% pull
    list_name$result[1, 8] <- list_name$rf_pred_2 %>% roc_auc(truth = class, paste0(".pred_", disease_name)) %>% dplyr::select(.estimate) %>% pull
  }
  
  return(list_name)
}


# ==============================================================================
# IV. Calibration Curve
# ==============================================================================

#' Generate Data for Calibration Curve
#'
#' This function calculates the Brier score and generates data for plotting a
#' calibration curve from prediction results.
#'
#' @param prediction_res A data frame containing prediction results, including true 
#'        class labels and predicted probabilities.
#' @param disease_name The name of the positive class (disease).
#'
#' @return A list containing the Brier score (`brier_score`) and a data frame 
#'         with data for plotting the calibration curve (`curve_data`).
#'         
do_calib_curve <- function(prediction_res, disease_name) {
  output_list <- list()
  
  # Prepare data for calibration
  calib_data <- prediction_res %>%
    dplyr::select(class, paste0(".pred_", disease_name)) %>%
    dplyr::rename(prob = paste0(".pred_", disease_name)) %>%
    mutate(observed = as.numeric(class == disease_name))
  
  # Calculate Brier score
  output_list$brier_score <- mean((calib_data$observed - calib_data$prob)^2)
  
  # Generate data for calibration curve
  output_list$curve_data <- calib_data %>%
    mutate(prob_bin = cut(prob, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(prob_bin) %>%
    summarise(mean_prob = mean(prob), mean_observed = mean(observed), .groups = 'drop')
  
  return(output_list)
}