# ==============================================================================
# WBT-DC Pipeline for Disease Classification
# ==============================================================================
#
# This script implements the WBT-DC (Whole-Blood Transcriptome-Derived Classifier)
# pipeline for classifying various diseases using gene expression data.
# It involves data loading, preprocessing (resampling, GSVA),
# machine learning (Random Forest), and evaluation (ROC curve analysis).
#
# Dependencies: tidyverse, tidymodels, doParallel, GSVA, ggplot2, ROSE, pROC, vip
#
# ==============================================================================

# Load necessary packages
pacman::p_load(tidyverse, tidymodels, doParallel, GSVA, ggplot2, ROSE, pROC, vip)

# Source pipeline_function.R for helper functions (assuming it's in the same directory or accessible)
# This line is commented out as the agent assumes the functions are already loaded or accessible.
# If not, it should be uncommented and potentially adjusted for the correct path.
# source("pipeline_function.R")

# ==============================================================================
# I. Crohn's Disease (CD) Analysis
# ==============================================================================

message("Starting Crohn's Disease (CD) Analysis...")

# ------------------------------------------------------------------------------
# 1. Data Loading and Annotation Preparation - CD
# ------------------------------------------------------------------------------

# Load differential gene expression (DEG) results for CD
cd_1016_deg <- readRDS("github/data/DEGs/cd_1016_deg.rds")
# Load gene expression matrices for CD
cd_expr_matrix <- readRDS("github/data/genes_expression/cd_expr_list.rds")
# Load annotation data for CD
cd_annot <- readRDS("github/data/annotation/cd_annot.rds")
# Load TPM (Transcripts Per Million) lists for CD (currently unused in snippet, but loaded)
cd_tpm_list <- readRDS("github/data/genes_expression/cd_tpm_list.rds")

# Prepare annotation data for different CD cohorts
cd_annot$data_87 <- cd_annot$data_87 %>%
  dplyr::select("ID", "class") %>%
  mutate(class = factor(class, levels = c("CD", "Control")))

cd_annot$data_array_235 <- cd_annot$data_array_235 %>%
  dplyr::select("ID", "class") %>%
  mutate(class = factor(class, levels = c("CD", "Control")))

# ------------------------------------------------------------------------------
# 2. Gene List Generation for WBT-DC Pipeline - CD
# ------------------------------------------------------------------------------

# Generate a list of gene sets (features) for GSVA for CD.
# This simulates selecting subsets of significant up-regulated and down-regulated genes.
cd_gene_list <- list()
set.seed(1234) # for reproducibility of gene list selection
for (i in 1:100) {
  cd_gene_list[[paste0("up_feature_", i)]] <- cd_1016_deg$DEG_res$sig_up$gene_id[1:200][sample(140, sample(61:100, 1))]
}
for (i in 1:100) {
  cd_gene_list[[paste0("down_feature_", i)]] <- cd_1016_deg$DEG_res$sig_down$gene_id[1:200][sample(140, sample(61:100, 1))]
}

# ------------------------------------------------------------------------------
# 3. Oversampling and GSVA Score Calculation (WBT-DC) - CD
# ------------------------------------------------------------------------------

# Perform oversampling on the primary CD training dataset (data_1016)
resample_cd_data <- list()
resample_cd_data$data_1016 <- get_resample_data(
  gene_expression_file = cd_expr_matrix$data_1016,
  annot_file = cd_annot$data_1016,
  seed = 5
)

# Calculate GSVA scores for oversampled training data and two testing datasets
cd_gsva_score <- list()
cd_gsva_score$data_1016_over <- do_gsva(
  tpm_matrix = resample_cd_data$data_1016$expr_over,
  gsva_feature_list = cd_gene_list,
  annot_matrix = resample_cd_data$data_1016$annot_over,
  disease_name = "CD",
  control_name = "Control"
)

cd_gsva_score$data_87_over <- do_gsva(
  tpm_matrix = cd_expr_matrix$data_87,
  gsva_feature_list = cd_gene_list,
  annot_matrix = cd_annot$data_87,
  disease_name = "CD",
  control_name = "Control"
)

cd_gsva_score$data_235_over <- do_gsva(
  tpm_matrix = cd_expr_matrix$data_235,
  gsva_feature_list = cd_gene_list,
  annot_matrix = cd_annot$data_array_235,
  disease_name = "CD",
  control_name = "Control"
)

# ------------------------------------------------------------------------------
# 4. Random Forest Model Training (WBT-DC Pipeline) - CD
# ------------------------------------------------------------------------------

# Train Random Forest model using GSVA scores from oversampled data
cd_ml_res <- list()
cd_ml_res$rf_over <- do_ml_rf(
  training_data = cd_gsva_score$data_1016_over,
  testing_data_1 = cd_gsva_score$data_87_over,
  testing_data_2 = cd_gsva_score$data_235_over,
  preprocess_step = "NON", # No additional preprocessing needed as GSVA scores are already features
  disease_name = "CD",
  seed = 13
)

# ------------------------------------------------------------------------------
# 5. Conventional Method (Top DEG Genes) - CD
# ------------------------------------------------------------------------------

# Select top 100 up and 100 down-regulated genes as features for conventional method
cd_deg_top100 <- c(cd_1016_deg$DEG_res$sig_up$gene_id[1:100], cd_1016_deg$DEG_res$sig_down$gene_id[1:100])

# Prepare data for conventional method training and testing
cd_1016_expr_top100 <- resample_cd_data$data_1016$expr_over[rownames(resample_cd_data$data_1016$expr_over) %in% cd_deg_top100, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(resample_cd_data$data_1016$annot_over) %>%
  relocate(class, .after = ID) %>%
  mutate(class = factor(class))

cd_87_expr_top100 <- cd_expr_matrix$data_87[rownames(cd_expr_matrix$data_87) %in% cd_deg_top100, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(cd_annot$data_87) %>%
  relocate(class, .after = ID) %>%
  mutate(class = factor(class))

cd_235_expr_top100 <- cd_expr_matrix$data_235[rownames(cd_expr_matrix$data_235) %in% cd_deg_top100, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(cd_annot$data_array_235) %>%
  relocate(class, .after = ID) %>%
  mutate(class = factor(class))

# Train Random Forest model using conventional DEG features
cd_ml_res$rf_conventional <- do_ml_rf(
  training_data = cd_1016_expr_top100,
  testing_data_1 = cd_87_expr_top100,
  testing_data_2 = cd_235_expr_top100,
  preprocess_step = "NON",
  disease_name = "CD",
  seed = 123
)

# ------------------------------------------------------------------------------
# 6. ROC Curve Analysis and Plotting - CD
# ------------------------------------------------------------------------------

# Prepare actual and predicted values for ROC curve generation
cd_auc <- list()
cd_auc$actural_87 <- ifelse(cd_ml_res$rf_over$rf_pred_1$class == "CD", 1, 0)
cd_auc$actural_235 <- ifelse(cd_ml_res$rf_over$rf_pred_2$class == "CD", 1, 0)
cd_auc$pre_new_87 <- cd_ml_res$rf_over$rf_pred_1$.pred_CD # WBT-DC predictions
cd_auc$pre_old_87 <- cd_ml_res$rf_conventional$rf_pred_1$.pred_CD # Conventional predictions
cd_auc$pre_new_235 <- cd_ml_res$rf_over$rf_pred_2$.pred_CD
cd_auc$pre_old_235 <- cd_ml_res$rf_conventional$rf_pred_2$.pred_CD

## Plot ROC for CD cohort 87
pdf("github/results/roc_comparison_cd_85.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(cd_auc$actural_87, cd_auc$pre_new_87,
                    main = "Crohn's Disease", percent = TRUE, col = "#DB9D85", identity = FALSE)
rocobj2 <- lines.roc(cd_auc$actural_87, cd_auc$pre_old_87, percent = TRUE, col = "black")
abline(a = 100, b = -1, col = "grey", lty = 2)
text(30, 35, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#DB9D85", cex = 0.8)
text(30, 30, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black", cex = 0.8)
legend("bottomright", legend = c("WBT-DC", "benchmark method"), col = c("#DB9D85", "black"), lwd = 2)
dev.off()

## Plot ROC for CD cohort 235
pdf("github/results/roc_comparison_cd_235.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(cd_auc$actural_235, cd_auc$pre_new_235,
                    main = "Crohn's Disease", percent = TRUE, col = "#DB9D85")
rocobj2 <- lines.roc(cd_auc$actural_235, cd_auc$pre_old_235, percent = TRUE, col = "black")
abline(a = 100, b = -1, col = "grey", lty = 2)
text(30, 25, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#DB9D85", cex = 0.8)
text(30, 20, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black", cex = 0.8)
legend("bottomright", legend = c("WBT-DC", "benchmark method"), col = c("#DB9D85", "black"), lwd = 2)
dev.off()

message("Crohn's Disease (CD) Analysis Completed.")

# ==============================================================================
# II. Ulcerative Colitis (UC) Analysis
# ==============================================================================

message("Starting Ulcerative Colitis (UC) Analysis...")

# ------------------------------------------------------------------------------
# 1. Data Loading and Annotation Preparation - UC
# ------------------------------------------------------------------------------

# Load differential gene expression (DEG) results for UC
uc_1016_deg <- readRDS("github/data/DEGs/uc_1016_deg.rds")
# Load gene expression matrices for UC
uc_expr_matrix <- readRDS("github/data/genes_expression/uc_expr_list.rds")
# Load annotation data for UC
uc_annot <- readRDS("github/data/annotation/uc_annot.rds")
# Load TPM (Transcripts Per Million) lists for UC (currently unused in snippet, but loaded)
uc_tpm_list <- readRDS("github/data/genes_expression/uc_tpm_list.rds")

# Prepare annotation data for different UC cohorts
uc_annot$data_87 <- uc_annot$data_87 %>%
  dplyr::select("ID", "class") %>%
  mutate(class = factor(class, levels = c("UC", "Control")))

uc_annot$data_array_235 <- uc_annot$data_array_235 %>%
  dplyr::select("ID", "class") %>%
  mutate(class = factor(class, levels = c("UC", "Control")))

# ------------------------------------------------------------------------------
# 2. Gene List Generation for WBT-DC Pipeline - UC
# ------------------------------------------------------------------------------

# Generate a list of gene sets (features) for GSVA for UC.
uc_gene_list <- list()
set.seed(12345) # for reproducibility of gene list selection
for (i in 1:100) {
  uc_gene_list[[paste0("up_feature_", i)]] <- uc_1016_deg$DEG_res$sig_up$gene_id[1:100][sample(100, sample(30:50, 1))]
}
for (i in 1:100) {
  uc_gene_list[[paste0("down_feature_", i)]] <- uc_1016_deg$DEG_res$sig_down$gene_id[1:100][sample(100, sample(30:50, 1))]
}

# ------------------------------------------------------------------------------
# 3. Oversampling and GSVA Score Calculation (WBT-DC) - UC
# ------------------------------------------------------------------------------

# Perform oversampling on the primary UC training dataset (data_1016)
resample_uc_data <- list()
resample_uc_data$data_1016 <- get_resample_data(
  gene_expression_file = uc_expr_matrix$data_1016,
  annot_file = uc_annot$data_1016,
  seed = 5
)

# Calculate GSVA scores for oversampled training data and two testing datasets
uc_gsva_score <- list()
uc_gsva_score$data_1016_over <- do_gsva(
  tpm_matrix = resample_uc_data$data_1016$expr_over,
  gsva_feature_list = uc_gene_list,
  annot_matrix = resample_uc_data$data_1016$annot_over,
  disease_name = "UC",
  control_name = "Control"
)

uc_gsva_score$data_87_over <- do_gsva(
  tpm_matrix = uc_expr_matrix$data_87,
  gsva_feature_list = uc_gene_list,
  annot_matrix = uc_annot$data_87,
  disease_name = "UC",
  control_name = "Control"
)

uc_gsva_score$data_235_over <- do_gsva(
  tpm_matrix = uc_expr_matrix$data_235,
  gsva_feature_list = uc_gene_list,
  annot_matrix = uc_annot$data_array_235,
  disease_name = "UC",
  control_name = "Control"
)

# ------------------------------------------------------------------------------
# 4. Random Forest Model Training (WBT-DC Pipeline) - UC
# ------------------------------------------------------------------------------

# Train Random Forest model using GSVA scores from oversampled data
uc_ml_res <- list()
uc_ml_res$rf_over <- do_ml_rf(
  training_data = uc_gsva_score$data_1016_over,
  testing_data_1 = uc_gsva_score$data_87_over,
  testing_data_2 = uc_gsva_score$data_235_over,
  preprocess_step = "NON",
  disease_name = "UC",
  seed = 123
)

# ------------------------------------------------------------------------------
# 5. Conventional Method (Top DEG Genes) - UC
# ------------------------------------------------------------------------------

# Select top 100 up and 100 down-regulated genes as features for conventional method
uc_deg_top100 <- c(uc_1016_deg$DEG_res$sig_up$gene_id[1:100], uc_1016_deg$DEG_res$sig_down$gene_id[1:100])

# Prepare data for conventional method training and testing
uc_1016_expr_top100 <- resample_uc_data$data_1016$expr_over[rownames(resample_uc_data$data_1016$expr_over) %in% uc_deg_top100, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(resample_uc_data$data_1016$annot_over) %>%
  relocate(class, .after = ID) %>%
  mutate(class = factor(class))

uc_87_expr_top100 <- uc_expr_matrix$data_87[rownames(uc_expr_matrix$data_87) %in% uc_deg_top100, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(uc_annot$data_87) %>%
  relocate(class, .after = ID) %>%
  mutate(class = factor(class))

uc_235_expr_top100 <- uc_expr_matrix$data_235[rownames(uc_expr_matrix$data_235) %in% uc_deg_top100, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(uc_annot$data_array_235) %>%
  relocate(class, .after = ID) %>%
  mutate(class = factor(class))

# Train Random Forest model using conventional DEG features
uc_ml_res$rf_conventional <- do_ml_rf(
  training_data = uc_1016_expr_top100,
  testing_data_1 = uc_87_expr_top100,
  testing_data_2 = uc_235_expr_top100,
  preprocess_step = "NON",
  disease_name = "UC",
  seed = 12
)

# ------------------------------------------------------------------------------
# 6. ROC Curve Analysis and Plotting - UC
# ------------------------------------------------------------------------------

# Prepare actual and predicted values for ROC curve generation
uc_auc <- list()
uc_auc$actural_87 <- ifelse(uc_ml_res$rf_over$rf_pred_1$class == "UC", 1, 0)
uc_auc$actural_235 <- ifelse(uc_ml_res$rf_over$rf_pred_2$class == "UC", 1, 0)
uc_auc$pre_new_87 <- uc_ml_res$rf_over$rf_pred_1$.pred_UC
uc_auc$pre_old_87 <- uc_ml_res$rf_conventional$rf_pred_1$.pred_UC
uc_auc$pre_new_235 <- uc_ml_res$rf_over$rf_pred_2$.pred_UC
uc_auc$pre_old_235 <- uc_ml_res$rf_conventional$rf_pred_2$.pred_UC

## Plot ROC for UC cohort 87
pdf("github/results/roc_comparison_uc_85.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(uc_auc$actural_87, uc_auc$pre_new_87,
                    main = "Ulcerative Colitis ", percent = TRUE, col = "#9DB469", identity = FALSE)
rocobj2 <- lines.roc(uc_auc$actural_87, uc_auc$pre_old_87, percent = TRUE, col = "black")
abline(a = 100, b = -1, col = "grey", lty = 2)
text(30, 25, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#9DB469", cex = 0.8)
text(30, 20, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black", cex = 0.8)
legend("bottomright", legend = c("WBT-DC", "benchmark method"), col = c("#9DB469", "black"), lwd = 2)
dev.off()

## Plot ROC for UC cohort 235
pdf("github/results/roc_comparison_uc_235.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(uc_auc$actural_235, uc_auc$pre_new_235,
                    main = "Ulcerative Colitis", percent = TRUE, col = "#9DB469", identity = FALSE)
rocobj2 <- lines.roc(uc_auc$actural_235, uc_auc$pre_old_235, percent = TRUE, col = "black")
abline(a = 100, b = -1, col = "grey", lty = 2)
text(30, 25, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#9DB469", cex = 0.8)
text(30, 20, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black", cex = 0.8)
legend("bottomright", legend = c("WBT-DC", "benchmark method"), col = c("#9DB469", "black"), lwd = 2)
dev.off()

message("Ulcerative Colitis (UC) Analysis Completed.")

# ==============================================================================
# III. Amyotrophic Lateral Sclerosis (ALS) Analysis
# ==============================================================================

message("Starting Amyotrophic Lateral Sclerosis (ALS) Analysis...")

# ------------------------------------------------------------------------------
# 1. Data Loading and Annotation Preparation - ALS
# ------------------------------------------------------------------------------

# Load differential gene expression (DEG) results for ALS
als_741_deg <- readRDS("github/data/DEGs/als_741_DEGs.rds")
# Load gene expression matrices for ALS
als_expr_matrix <- readRDS("github/data/genes_expression/als_expr_list.rds")
# Load annotation data for ALS
als_annot <- readRDS("github/data/annotation/als_annot.rds")

# Prepare annotation data for ALS cohort
als_annot$data_85 <- als_annot$data_85 %>%
  dplyr::select("ID", "class") %>%
  mutate(class = factor(class)) %>%
  mutate(class = factor(class, levels = c("ALS", "Control")))

# ------------------------------------------------------------------------------
# 2. Gene List Generation for WBT-DC Pipeline - ALS
# ------------------------------------------------------------------------------

# Generate a list of gene sets (features) for GSVA for ALS.
als_gene_list <- list()
set.seed(1235) # for reproducibility of gene list selection
for (i in 1:100) {
  als_gene_list[[paste0("up_feature_", i)]] <- als_741_deg$up_gene[1:160][sample(160, sample(51:100, 1))]
}
for (i in 1:100) {
  als_gene_list[[paste0("down_feature_", i)]] <- als_741_deg$down_gene[1:160][sample(160, sample(51:100, 1))]
}

# ------------------------------------------------------------------------------
# 3. Oversampling and GSVA Score Calculation (WBT-DC) - ALS
# ------------------------------------------------------------------------------

# Perform oversampling on the primary ALS training dataset (data_741)
resample_als_data <- list()
resample_als_data$data_741 <- get_resample_data(
  gene_expression_file = als_expr_matrix$data_741,
  annot_file = als_annot$data_741,
  seed = 5
)

# Calculate GSVA scores for oversampled training data and one testing dataset
als_gsva_score <- list()
als_gsva_score$data_741_over <- do_gsva(
  tpm_matrix = resample_als_data$data_741$expr_over,
  gsva_feature_list = als_gene_list,
  annot_matrix = resample_als_data$data_741$annot_over,
  disease_name = "ALS",
  control_name = "Control"
)

als_gsva_score$data_85_over <- do_gsva(
  tpm_matrix = als_expr_matrix$data_85,
  gsva_feature_list = als_gene_list,
  annot_matrix = als_annot$data_85,
  disease_name = "ALS",
  control_name = "Control"
)

# ------------------------------------------------------------------------------
# 4. Random Forest Model Training (WBT-DC Pipeline) - ALS
# ------------------------------------------------------------------------------

# Train Random Forest model using GSVA scores from oversampled data
als_ml_res <- list()
als_ml_res$rf_over <- do_ml_rf(
  training_data = als_gsva_score$data_741_over,
  testing_data_1 = als_gsva_score$data_85_over,
  testing_data_2 = NULL, # Only one testing dataset for ALS
  preprocess_step = "NON",
  disease_name = "ALS",
  seed = 123
)

# ------------------------------------------------------------------------------
# 5. Conventional Method (Top DEG Genes) - ALS
# ------------------------------------------------------------------------------

# Select top 100 up and 100 down-regulated genes as features for conventional method
als_deg_top100 <- c(als_741_deg$up_gene[1:100], als_741_deg$down_gene[1:100])

# Prepare data for conventional method training and testing
als_741_expr_top100 <- resample_als_data$data_741$expr_over[rownames(resample_als_data$data_741$expr_over) %in% als_deg_top100, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(resample_als_data$data_741$annot_over %>% dplyr::select(ID, class)) %>%
  relocate(class, .after = ID) %>%
  mutate(class = factor(class, levels = c("ALS", "Control")))

als_85_expr_top100 <- als_expr_matrix$data_85[rownames(als_expr_matrix$data_85) %in% als_deg_top100, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  left_join(als_annot$data_85 %>% dplyr::select(ID, class)) %>%
  relocate(class, .after = ID) %>%
  mutate(class = factor(class, levels = c("ALS", "Control")))

# Train Random Forest model using conventional DEG features
als_ml_res$rf_convention <- do_ml_rf(
  training_data = als_741_expr_top100,
  testing_data_1 = als_85_expr_top100,
  testing_data_2 = NULL,
  preprocess_step = "NON",
  disease_name = "ALS",
  seed = 123
)

# ------------------------------------------------------------------------------
# 6. ROC Curve Analysis and Plotting - ALS
# ------------------------------------------------------------------------------

# Prepare actual and predicted values for ROC curve generation
als_auc <- list()
als_auc$actural_85_1 <- ifelse(als_ml_res$rf_over$rf_pred_1$class == "ALS", 1, 0)
als_auc$pre_new_85 <- als_ml_res$rf_over$rf_pred_1$.pred_ALS # WBT-DC predictions
als_auc$actural_85_2 <- ifelse(als_ml_res$rf_convention$rf_pred_1$class == "ALS", 1, 0)
als_auc$pre_old_85 <- als_ml_res$rf_convention$rf_pred_1$.pred_ALS # Conventional predictions

## Plot ROC for ALS cohort 85
pdf("github/results/roc_comparison_als_85.pdf", width = 6, height = 6)
rocobj1 <- plot.roc(als_auc$actural_85_1, als_auc$pre_new_85,
                    main = "Amyotrophic Lateral Sclerosis ", percent = TRUE, col = "#BB9FE0")
rocobj2 <- lines.roc(als_auc$actural_85_2, als_auc$pre_old_85, percent = TRUE, col = "black")
abline(a = 100, b = -1, col = "grey", lty = 2)
text(30, 25, labels = paste("AUC =", round(auc(rocobj1), 2)), col = "#BB9FE0", cex = 0.8)
text(30, 20, labels = paste("AUC =", round(auc(rocobj2), 2)), col = "black", cex = 0.8)
legend("bottomright", legend = c("WBT-DC", "benchmark method"), col = c("#BB9FE0", "black"), lwd = 2)
dev.off()

message("Amyotrophic Lateral Sclerosis (ALS) Analysis Completed.")

# ==============================================================================
# IV. Rheumatoid Arthritis (RA) Analysis
# ==============================================================================

message("Starting Rheumatoid Arthritis (RA) Analysis...")

# ------------------------------------------------------------------------------
# 1. Data Loading and Annotation Preparation - RA
# ------------------------------------------------------------------------------

# Load differential gene expression (DEG) results for RA
ra_sza_deg <- readRDS("github/data/DEGs/ra_sza_deg.rds")
# Load gene expression matrices for RA
ra_expr_matrix <- readRDS("github/data/genes_expression/ra_expr_list.rds")
# Load annotation data for RA
ra_annot <- readRDS("github/data/annotation/ra_annot.rds")

# Prepare annotation data for RA cohort
ra_annot$data_275 <- ra_annot$data_275 %>%
  dplyr::select("ID", "class") %>%
  mutate(class = factor(class)) %>%
  mutate(class = factor(class, levels = c("RA", "Control")))

# ------------------------------------------------------------------------------
# 2. Gene List Generation for WBT-DC Pipeline - RA
# ------------------------------------------------------------------------------

# Generate a list of gene sets (features) for GSVA for RA.
ra_gene_list <- list()
set.seed(1234) # for reproducibility of gene list selection
for (i in 1:100) {
  ra_gene_list[[paste0("up_feature_", i)]] <- ra_sza_deg$DEG_res$sig_up$gene_id[1:200][sample(200, sample(61:100, 1))]
}
for (i in 1:100) {
  ra_gene_list[[paste0("down_feature_", i)]] <- ra_sza_deg$DEG_res$sig_down$gene_id[1:200][sample(200, sample(61:100, 1))]
}

# ------------------------------------------------------------------------------
# 3. Oversampling and GSVA Score Calculation (WBT-DC) - RA
# ------------------------------------------------------------------------------

# Perform oversampling on the primary RA training dataset (sza)
resample_ra_data <- list()
resample_ra_data$sza <- get_resample_data(
  gene_expression_file = ra_expr_matrix$sza,
  annot_file = ra_annot$sza,
  seed = 5
)

# Calculate GSVA scores for oversampled training data and one testing dataset
ra_gsva_score <- list()
ra_gsva_score$data_sza_over <- do_gsva(
  tpm_matrix = resample_ra_data$sza$expr_over,
  gsva_feature_list = ra_gene_list,
  annot_matrix = resample_ra_data$sza$annot_over,
  disease_name = "RA",
  control_name = "Control"
)

ra_gsva_score$data_275_over <- do_gsva(
  tpm_matrix = ra_expr_matrix$data_275,
  gsva_feature_list = ra_gene_list,
  annot_matrix = ra_annot$data_275,
  disease_name = "RA",
  control_name = "Control"
)

# ------------------------------------------------------------------------------
# 4. Random Forest Model Training (WBT-DC Pipeline) - RA
# ------------------------------------------------------------------------------

# Train Random Forest model using GSVA scores from oversampled data
ra_ml_res <- list()
ra_ml_res$rf_over <- do_ml_rf(
  training_data = ra_gsva_score$data_sza_over,
  testing_data_1 = ra_gsva_score$data_275_over,
  testing_data_2 = NULL, # Only one testing dataset for RA
  preprocess_step = "NON",
  disease_name = "RA",
  seed = 13
)

# ------------------------------------------------------------------------------
# 5. ROC Curve Analysis and Plotting - RA
# ------------------------------------------------------------------------------

# Prepare actual and predicted values for ROC curve generation
ra_auc <- list()
ra_auc$actural_275 <- ifelse(ra_ml_res$rf_over$rf_pred_1$class == "RA", 1, 0)
ra_auc$pre_275 <- ra_ml_res$rf_over$rf_pred_1$.pred_RA # WBT-DC predictions

## Plot ROC for RA cohort 275
pdf("github/results/roc_ra_275.pdf", width = 6, height = 6)
rocobj <- plot.roc(ra_auc$actural_275, ra_auc$pre_275,
                   main = "Confidence intervals",
                   percent = TRUE,
                   ci = TRUE,
                   print.auc = TRUE)
ciobj <- ci.se(rocobj,
               specificities = seq(0, 100, 5))
plot(ciobj, type = "shape", col = "#1c61b6AA")
plot(ci(rocobj, of = "thresholds", thresholds = "best"))
dev.off()

message("Rheumatoid Arthritis (RA) Analysis Completed.")