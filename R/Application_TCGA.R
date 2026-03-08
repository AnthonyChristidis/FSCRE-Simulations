# ---------------------------------------------------------------------------
# Application 2: Proteogenomics (TCGA BRCA) - Protein Abundance Prediction
# ---------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# --- 1. Load Required Libraries ---
library(curatedTCGAData)
library(TCGAutils)
library(glmnet)
library(cellWise)
library(randomGLM)
library(caret)
library(srlars) 
library(regcell)

# Load competitor scripts (adjust paths as needed)
source("R/SparseShootingS/sparseShootingS.R")

set.seed(456)

# --- 2. Data Downloading and Processing ---
cat("\n--- Preparing TCGA BRCA Dataset ---\n")

# Create data directory if it doesn't exist
if (!dir.exists("data")) dir.create("data")

# We want RNA-seq (gene expression) and RPPA (Reverse Phase Protein Array) for Breast Cancer
# RNA: "BRCA_RNASeq2GeneNorm-20160128"
# Protein: "BRCA_RPPAArray-20160128"

# Path to save the processed matched data
data_file <- "data/TCGA_BRCA_Matched.rds"

if(!file.exists(data_file)) {
  cat("Downloading TCGA BRCA data (may take a minute)...\n")
  
  # Download the specific assays
  # Note: curatedTCGAData manages its own internal cache (usually in ExperimentHub)
  # We just extract the specific matrices we need and save those to our local data/ folder
  brca_data <- curatedTCGAData(
    diseaseCode = "BRCA",
    assays = c("RNASeq2GeneNorm", "RPPAArray"),
    dry.run = FALSE
  )
  
  # TCGAutils helps match samples that have BOTH RNA and Protein data
  matched_data <- intersectOperations(brca_data)
  
  # Extract the matrices
  # Transpose so rows are patients, columns are genes/proteins
  rna_mat <- t(assay(matched_data[["BRCA_RNASeq2GeneNorm-20160128"]]))
  prot_mat <- t(assay(matched_data[["BRCA_RPPAArray-20160128"]]))
  
  # Clean up patient IDs so they match exactly
  # TCGA barcodes have format "TCGA-XX-XXXX-01A-..." we just need the patient ID part
  rownames(rna_mat) <- substr(rownames(rna_mat), 1, 12)
  rownames(prot_mat) <- substr(rownames(prot_mat), 1, 12)
  
  # Find exact intersection again just to be safe
  common_patients <- intersect(rownames(rna_mat), rownames(prot_mat))
  
  rna_mat <- rna_mat[common_patients, ]
  prot_mat <- prot_mat[common_patients, ]
  
  # Save the finalized, matched matrices to avoid re-downloading and re-processing
  saveRDS(list(rna = rna_mat, prot = prot_mat), data_file)
  
} else {
  cat("Loading TCGA data from local file...\n")
  mats <- readRDS(data_file)
  rna_mat <- mats$rna
  prot_mat <- mats$prot
}

# --- 3. Define the Prediction Task ---

# Response (y): Protein abundance of a key driver.
# ESR1 (Estrogen Receptor) is highly relevant in Breast Cancer.
target_protein <- "ER.alpha" # This is the RPPA name for ESR1

if (!(target_protein %in% colnames(prot_mat))) {
  stop("Target protein not found. Available proteins: ", paste(head(colnames(prot_mat)), collapse=", "))
}

y_full <- prot_mat[, target_protein]

# Remove any NAs in the target protein
valid_y_idx <- which(!is.na(y_full))
y <- y_full[valid_y_idx]
X <- rna_mat[valid_y_idx, ]

# Filter Predictors (X):
# 1. Remove genes with zero variance or very low expression
gene_vars <- apply(X, 2, var)
valid_genes <- which(gene_vars > 1e-4)
X <- X[, valid_genes]

# 2. Filter to top 1000 most variable genes to create p >> n setting
gene_vars <- apply(X, 2, var)
top_genes <- order(gene_vars, decreasing = TRUE)[1:1000]
X <- X[, top_genes]

# Standardize predictors
X <- scale(X)

cat(sprintf("Final Data Dimensions: n = %d, p = %d\n", nrow(X), ncol(X)))

# --- 4. Define Model Evaluation Helper ---
# (This is identical to Application 1)
run_models <- function(x_tr, y_tr, x_te, y_te, prefix) {
  
  p <- ncol(x_tr)
  mspe_results <- c()
  
  # 1. Elastic Net (EN)
  fit_en <- tryCatch({
    fit <- cv.glmnet(x = x_tr, y = y_tr, alpha = 0.5)
    mean((predict(fit, x_te, s = "lambda.min") - y_te)^2)
  }, error = function(e) NA)
  mspe_results["ElasticNet"] <- fit_en
  
  # Shared DDC Imputation for Baselines
  x_imp <- x_tr
  y_imp <- y_tr
  tryCatch({
    ddc_out <- cellWise::DDC(cbind(x_tr, y_tr), DDCpars = list(fastDDC = TRUE, silent = TRUE))
    x_imp <- ddc_out$Ximp[, 1:p]
    y_imp <- ddc_out$Ximp[, p+1]
  }, error = function(e) NULL)
  
  # 2. DDC + EN
  fit_ddc_en <- tryCatch({
    fit <- cv.glmnet(x = x_imp, y = y_imp, alpha = 0.5)
    mean((predict(fit, x_te, s = "lambda.min") - y_te)^2)
  }, error = function(e) NA)
  mspe_results["DDC_EN"] <- fit_ddc_en
  
  # 3. DDC + Random GLM
  fit_ddc_rglm <- tryCatch({
    fit <- randomGLM(x_imp, y_imp, classify = FALSE, nBags = 100, keepModels = TRUE, nThreads = 1, verbose = 0)
    mean((predict(fit, newdata = x_te) - y_te)^2)
  }, error = function(e) NA)
  mspe_results["DDC_RGLM"] <- fit_ddc_rglm
  
  # 4. Sparse Shooting S
  fit_sps <- tryCatch({
    fit <- sparseshooting(x = x_tr, y = y_tr, wvalue = 3, nlambda = 50)
    preds <- fit$coef[1] + x_te %*% fit$coef[-1]
    mean((preds - y_te)^2)
  }, error = function(e) NA)
  mspe_results["Sparse_S"] <- fit_sps
  
  # 5. CR-Lasso
  fit_crlasso <- tryCatch({
    fit <- regcell::sregcell_std(y_tr, x_tr)
    preds <- fit$intercept_hat + x_te %*% fit$betahat
    mean((preds - y_te)^2)
  }, error = function(e) NA)
  mspe_results["CR_Lasso"] <- fit_crlasso
  
  # 6. RLARS (Proposed, K=1)
  fit_rlars <- tryCatch({
    fit <- srlars(x_tr, y_tr, n_models = 1, tolerance = 0.01, robust = TRUE, compute_coef = TRUE)
    mean((predict(fit, newx = x_te, dynamic = FALSE) - y_te)^2)
  }, error = function(e) NA)
  mspe_results["RLARS"] <- fit_rlars
  
  # 7. FSCRE (Proposed, K=10)
  fit_fscre <- tryCatch({
    fit <- srlars(x_tr, y_tr, n_models = 10, tolerance = 0.01, robust = TRUE, compute_coef = TRUE)
    mean((predict(fit, newx = x_te, dynamic = FALSE) - y_te)^2)
  }, error = function(e) NA)
  mspe_results["FSCRE"] <- fit_fscre
  
  # Prepend prefix to names
  names(mspe_results) <- paste0(prefix, "_", names(mspe_results))
  return(mspe_results)
}

# --- 5. Simulation Loop ---
N_splits <- 50
results_list <- list()

cat("\n--- Starting TCGA Data Splits ---\n")

for (i in 1:N_splits) {
  cat(sprintf("Split %d / %d...\n", i, N_splits))
  
  # 70/30 Train/Test Split
  train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
  
  x_train <- X[train_idx, ]
  y_train <- y[train_idx]
  x_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # A. Run on ORIGINAL Data
  mspe_orig <- run_models(x_train, y_train, x_test, y_test, prefix = "Orig")
  
  # B. Introduce Artificial Contamination to Training Data (10% marginal shift)
  # Using 10% here to show it can handle heavier contamination than the GDSC example
  x_train_cont <- x_train
  n_cells <- nrow(x_train) * ncol(x_train)
  contam_idx <- sample(1:n_cells, round(0.10 * n_cells))
  x_train_cont[contam_idx] <- runif(length(contam_idx), min = 5, max = 15) 
  
  # Run on CONTAMINATED Data
  mspe_cont <- run_models(x_train_cont, y_train, x_test, y_test, prefix = "Contam")
  
  # Store results
  results_list[[i]] <- c(mspe_orig, mspe_cont)
}

# --- 6. Summarize and Save Results ---
final_results_df <- do.call(rbind, results_list)

cat("\n--- FINAL TCGA RESULTS (Average MSPE) ---\n")
print(round(colMeans(final_results_df, na.rm = TRUE), 4))

if (!dir.exists("results")) dir.create("results")
saveRDS(final_results_df, "results/Application_TCGA_Results.rds")
cat("\nResults saved to results/Application2_TCGA_Results.rds\n")