# ------------------------------------------------------------------------
# Application: Proteogenomics (TCGA BRCA) - Protein Abundance Prediction
# ------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(curatedTCGAData)
library(TCGAutils)
library(glmnet)
library(cellWise)
library(randomGLM)
library(caret)
library(srlars) 
library(regcell)

# Load required source files
source("SparseShootingS/sparseShootingS.R")

# ------------------------------------------------------------------------

# Set seed
set.seed(0)

# _________________________________
# Data Downloading and Processing 
# _________________________________

cat("\n--- Preparing TCGA BRCA Dataset ---\n")

if (!dir.exists("data")) dir.create("data")
if (!dir.exists("results")) dir.create("results")

data_file <- "data/TCGA_BRCA_Matched.rds"

if(!file.exists(data_file)) {
  cat("Downloading TCGA BRCA data (may take a minute)...\n")
  
  # Download the specific assays. 
  # version="2.1.1" is the latest stable release
  brca_data <- curatedTCGAData::curatedTCGAData(
    diseaseCode = "BRCA",
    assays = c("RNASeq2GeneNorm", "RPPAArray"),
    version = "2.1.1",
    dry.run = FALSE
  )
  
  # TCGAutils helps match samples that have BOTH RNA and Protein data
  matched_data <- MultiAssayExperiment::intersectColumns(brca_data)
  
  # Extract the matrices. 
  assay_names <- names(experiments(matched_data))
  rna_name <- assay_names[grep("RNASeq2GeneNorm", assay_names)]
  prot_name <- assay_names[grep("RPPAArray", assay_names)]
  
  rna_mat <- t(assay(matched_data[[rna_name]]))
  prot_mat <- t(assay(matched_data[[prot_name]]))
  
  # Clean up patient IDs so they match exactly
  rownames(rna_mat) <- substr(rownames(rna_mat), 1, 12)
  rownames(prot_mat) <- substr(rownames(prot_mat), 1, 12)
  
  # Handle potential replicates
  rna_mat <- rna_mat[!duplicated(rownames(rna_mat)), ]
  prot_mat <- prot_mat[!duplicated(rownames(prot_mat)), ]
  
  # Find exact intersection
  common_patients <- intersect(rownames(rna_mat), rownames(prot_mat))
  
  rna_mat <- rna_mat[common_patients, ]
  prot_mat <- prot_mat[common_patients, ]
  
  cat(sprintf("Matched %d patients with both RNA and Protein data.\n", length(common_patients)))
  
  saveRDS(list(rna = rna_mat, prot = prot_mat), data_file)
  
} else {
  cat("Loading TCGA data from local file...\n")
  mats <- readRDS(data_file)
  rna_mat <- mats$rna
  prot_mat <- mats$prot
}


# _______________________________
# Define Model Evaluation Helper 
# _______________________________

runModels <- function(x_tr, y_tr, x_te, y_te, prefix) {
  
  p <- ncol(x_tr)
  mspe_results <- c()
  
  # 1. Elastic Net (EN)
  fit_en <- tryCatch({
    fit <- cv.glmnet(x = x_tr, y = y_tr, alpha = 3/4)
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
    fit <- cv.glmnet(x = x_imp, y = y_imp, alpha = 3/4)
    mean((predict(fit, x_te, s = "lambda.min") - y_te)^2)
  }, error = function(e) NA)
  mspe_results["DDC_EN"] <- fit_ddc_en
  
  # 3. DDC + Random GLM
  fit_ddc_rglm <- tryCatch({
    fit <- randomGLM(x_imp, y_imp, classify = FALSE, nBags = 100, keepModels = TRUE, nThreads = 1, verbose = 0)
    mean((predict(fit, newdata = x_te) - y_te)^2)
  }, error = function(e) NA)
  mspe_results["DDC_RGLM"] <- fit_ddc_rglm
  
  # # 4. Sparse Shooting S
  # fit_sps <- tryCatch({
  #   fit <- sparseshooting(x = x_tr, y = y_tr, wvalue = 3, nlambda = 50)
  #   preds <- fit$coef[1] + x_te %*% fit$coef[-1]
  #   mean((preds - y_te)^2)
  # }, error = function(e) NA)
  # mspe_results["Sparse_S"] <- fit_sps
  
  # # 5. CR-Lasso
  # fit_crlasso <- tryCatch({
  #   fit <- regcell::sregcell_std(y_tr, x_tr)
  #   preds <- fit$intercept_hat + x_te %*% fit$betahat
  #   mean((preds - y_te)^2)
  # }, error = function(e) NA)
  # mspe_results["CR_Lasso"] <- fit_crlasso
  
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


# ________________________________
# OUTER LOOP OVER TARGET PROTEINS
# ________________________________

proteins_to_run <- c("ER-alpha", "EGFR")
N_splits <- 50

for (target_protein in proteins_to_run) {
  
  cat(sprintf("\n\n*******************************************************\n"))
  cat(sprintf("STARTING ANALYSIS FOR TARGET PROTEIN: %s\n", target_protein))
  cat(sprintf("*******************************************************\n"))
  
  # _______________________
  # Define Prediction Task 
  # _______________________
  
  if (!(target_protein %in% colnames(prot_mat))) {
    warning("Target protein ", target_protein, " not found. Skipping.")
    next
  }
  
  y_full <- prot_mat[, target_protein]
  
  # Remove any NAs in the target protein
  valid_y_idx <- which(!is.na(y_full))
  y <- y_full[valid_y_idx]
  X_full <- rna_mat[valid_y_idx, ]
  
  # Filter Predictors (X):
  # 1. Remove genes with zero variance or very low expression
  gene_vars <- apply(X_full, 2, var)
  valid_genes <- which(gene_vars > 1e-4)
  X_filtered <- X_full[, valid_genes]
  
  # 2. Correlation filter to ensure strong baseline signal
  cat("Filtering to top 500 most correlated genes...\n")
  gene_corrs <- abs(apply(X_filtered, 2, function(x) cor(x, y)))
  gene_corrs[is.na(gene_corrs)] <- 0
  top_genes <- order(gene_corrs, decreasing = TRUE)[1:500]
  X <- X_filtered[, top_genes]
  
  # Standardize predictors and response
  X <- scale(X)
  y <- scale(y) 
  
  cat(sprintf("Data Dimensions for %s: n = %d, p = %d\n", target_protein, nrow(X), ncol(X)))
  
  # ________________
  # Simulation Loop
  # ________________
  
  results_list <- list()
  
  for (i in 1:N_splits) {
    cat(sprintf("  Split %d / %d...\n", i, N_splits))
    
    # Train/Test Split (Note: p=0.1 means 10% in test, 90% in train. Adjust if you meant 70/30)
    train_idx <- createDataPartition(y, p = 0.7, list = FALSE) # Changed to 0.7 for standard 70/30 split
    
    x_train <- X[train_idx, ]
    y_train <- y[train_idx]
    x_test <- X[-train_idx, ]
    y_test <- y[-train_idx]
    
    # A. Run on ORIGINAL Data
    mspe_orig <- runModels(x_train, y_train, x_test, y_test, prefix = "Orig")
    
    # B. Introduce TARGETED Artificial Contamination 
    x_train_cont <- x_train
    
    # 1. Identify the "most important" variables from a quick, clean EN fit
    quick_en <- cv.glmnet(x_train, y_train, alpha = 3/4)
    clean_coefs <- as.numeric(coef(quick_en, s = "lambda.min"))[-1]
    
    # Get the indices of the top 10 most important genes
    n_top_genes <- min(10, sum(clean_coefs != 0))
    if (n_top_genes == 0) n_top_genes <- 10 
    top_var_indices <- order(abs(clean_coefs), decreasing = TRUE)[1:n_top_genes]
    
    # 2. Poison ONLY these important variables
    # Corrupt 15% of the cells within these specific important columns
    n_to_corrupt <- round(0.15 * nrow(x_train))
    
    for (j in top_var_indices) {
        corrupt_rows <- sample(1:nrow(x_train), n_to_corrupt)
        # Data is scaled (mean 0, SD 1), so +/- 10 is a massive outlier
        x_train_cont[corrupt_rows, j] <- sample(c(10, -10), n_to_corrupt, replace = TRUE)
    }
    
    # Run on CONTAMINATED Data
    mspe_cont <- runModels(x_train_cont, y_train, x_test, y_test, prefix = "Contam")
    
    # Store results
    results_list[[i]] <- c(mspe_orig, mspe_cont)
  }
  
  # ___________________________
  # Summarize and Save Results
  # ___________________________
  
  final_results_df <- do.call(rbind, results_list)
  
  cat(sprintf("\n--- FINAL RESULTS for %s (Average MSPE) ---\n", target_protein))
  print(round(colMeans(final_results_df, na.rm = TRUE), 4))
  
  # Clean filename (replace hyphens or spaces to avoid saving issues)
  safe_protein_name <- gsub("-", "_", target_protein)
  filename <- sprintf("results/Application2_TCGA_%s_Results.rds", safe_protein_name)
  
  saveRDS(final_results_df, filename)
  cat(sprintf("Results saved to %s\n", filename))

  # ___________________________________________________________
  # Final Model Fit on Full Data (for Gene Selection Analysis)
  # ___________________________________________________________

  cat(sprintf("\n--- Fitting Final Models on Full Data for %s ---\n", target_protein))

  # A. Fit on Full Original Data
  fit_orig_fscre <- tryCatch({
      srlars(X, y, n_models = 10, tolerance = 0.01, robust = TRUE)
  }, error = function(e) list(active.sets = list()))
  
  orig_selected_genes <- colnames(X)[unique(unlist(fit_orig_fscre$active.sets))]

  # B. Create TARGETED Contaminated Version of the Full Data
  X_cont_full <- X
  
  # Identify important variables on full clean data
  quick_en_full <- cv.glmnet(X, y, alpha = 3/4)
  clean_coefs_full <- as.numeric(coef(quick_en_full, s = "lambda.min"))[-1]
  n_top_genes_full <- min(10, sum(clean_coefs_full != 0))
  if (n_top_genes_full == 0) n_top_genes_full <- 10
  top_var_indices_full <- order(abs(clean_coefs_full), decreasing = TRUE)[1:n_top_genes_full]
  
  # Poison those variables
  n_to_corrupt_full <- round(0.15 * nrow(X))
  for (j in top_var_indices_full) {
      corrupt_rows <- sample(1:nrow(X), n_to_corrupt_full)
      X_cont_full[corrupt_rows, j] <- sample(c(10, -10), n_to_corrupt_full, replace = TRUE)
  }

  # C. Fit Models on Contaminated Data
  fit_cont_fscre <- tryCatch({
      srlars(X_cont_full, y, n_models = 10, tolerance = 0.01, robust = TRUE)
  }, error = function(e) list(active.sets = list()))
  
  cont_selected_genes_fscre <- colnames(X)[unique(unlist(fit_cont_fscre$active.sets))]

  fit_cont_en <- tryCatch({
      cv.glmnet(x = X_cont_full, y = y, alpha = 3/4)
  }, error = function(e) NULL)
  
  if (!is.null(fit_cont_en)) {
      cont_coefs_en <- as.numeric(coef(fit_cont_en, s = "lambda.min"))[-1]
      cont_selected_genes_en <- colnames(X)[which(cont_coefs_en != 0)]
  } else {
      cont_selected_genes_en <- character(0)
  }

  # D. Save the Lists for Literature Checking
  genes_list <- list(
      FSCRE_Orig = orig_selected_genes, 
      FSCRE_Cont = cont_selected_genes_fscre, 
      EN_Cont = cont_selected_genes_en
  )
  
  genes_filename <- sprintf("results/Application2_TCGA_%s_Genes.rds", safe_protein_name)
  saveRDS(genes_list, genes_filename)
  
  cat(sprintf("Gene selection lists saved to %s\n", genes_filename))
}

