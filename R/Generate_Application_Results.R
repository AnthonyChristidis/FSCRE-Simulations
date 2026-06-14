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
library(R.utils)

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
  
  brca_data <- curatedTCGAData::curatedTCGAData(
    diseaseCode = "BRCA",
    assays = c("RNASeq2GeneNorm", "RPPAArray"),
    version = "2.1.1",
    dry.run = FALSE
  )
  
  matched_data <- MultiAssayExperiment::intersectColumns(brca_data)
  
  assay_names <- names(experiments(matched_data))
  rna_name <- assay_names[grep("RNASeq2GeneNorm", assay_names)]
  prot_name <- assay_names[grep("RPPAArray", assay_names)]
  
  rna_mat <- t(assay(matched_data[[rna_name]]))
  prot_mat <- t(assay(matched_data[[prot_name]]))
  
  rownames(rna_mat) <- substr(rownames(rna_mat), 1, 12)
  rownames(prot_mat) <- substr(rownames(prot_mat), 1, 12)
  
  rna_mat <- rna_mat[!duplicated(rownames(rna_mat)), ]
  prot_mat <- prot_mat[!duplicated(rownames(prot_mat)), ]
  
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
  selected_vars <- list()
  
  # 1. Elastic Net (EN)
  tryCatch({
    fit <- cv.glmnet(x = x_tr, y = y_tr, alpha = 3/4)
    mspe_results["ElasticNet"] <- mean((predict(fit, x_te, s = "lambda.min") - y_te)^2)
    coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
    selected_vars[["ElasticNet"]] <- colnames(x_tr)[which(coefs != 0)]
  }, error = function(e) {
    mspe_results["ElasticNet"] <<- NA
    selected_vars[["ElasticNet"]] <<- character(0)
  })
  
  # --- Shared Data Cleaning for DDC baselines (X-Only, No Leakage) ---
  x_imp <- x_tr
  y_imp <- y_tr
  x_te_imp <- x_te
  
  tryCatch({
    # Clean predictors ONLY
    ddc_out <- cellWise::DDC(x_tr, DDCpars = list(fastDDC = TRUE, silent = TRUE))
    x_imp <<- ddc_out$Ximp
    
    # Apply DDC model to test set
    ddc_test <- cellWise::DDCpredict(x_te, ddc_out)
    x_te_imp <<- ddc_test$Ximp
    
    # Univariate wrap for response
    y_imp <<- as.numeric(cellWise::wrap(as.matrix(y_tr))$Xw)
  }, error = function(e) {})
  
  # 2. DDC + EN
  tryCatch({
    fit <- cv.glmnet(x = x_imp, y = y_imp, alpha = 3/4)
    mspe_results["DDC_EN"] <- mean((predict(fit, x_te_imp, s = "lambda.min") - y_te)^2)
    coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
    selected_vars[["DDC_EN"]] <- colnames(x_imp)[which(coefs != 0)]
  }, error = function(e) {
    mspe_results["DDC_EN"] <<- NA
    selected_vars[["DDC_EN"]] <<- character(0)
  })
  
  # 3. DDC + Random GLM
  tryCatch({
    fit <- randomGLM(x_imp, y_imp, classify = FALSE, nBags = 100, keepModels = TRUE, nThreads = 1, verbose = 0)
    mspe_results["DDC_RGLM"] <- mean((predict(fit, newdata = x_te_imp) - y_te)^2)
    sel_counts <- colSums(fit$timesSelectedByForwardRegression)
    selected_vars[["DDC_RGLM"]] <- colnames(x_imp)[which(sel_counts > 0)]
  }, error = function(e) {
    mspe_results["DDC_RGLM"] <<- NA
    selected_vars[["DDC_RGLM"]] <<- character(0)
  })
  
  # 4. Sparse Shooting S
  tryCatch({
    fit <- sparseshooting(x = x_tr, y = y_tr, wvalue = 3, nlambda = 50)
    preds <- fit$coef[1] + x_te %*% fit$coef[-1]
    mspe_results["Sparse_S"] <- mean((preds - y_te)^2)
    coefs <- fit$coef[-1]
    selected_vars[["Sparse_S"]] <- colnames(x_tr)[which(coefs != 0)]
  }, error = function(e) {
    mspe_results["Sparse_S"] <<- NA
    selected_vars[["Sparse_S"]] <<- character(0)
  })
  
  # 5. CR-Lasso
  tryCatch({
    fit <- regcell::sregcell_std(y_tr, x_tr)
    preds <- fit$intercept_hat + x_te %*% fit$betahat
    mspe_results["CR_Lasso"] <- mean((preds - y_te)^2)
    coefs <- fit$betahat
    selected_vars[["CR_Lasso"]] <- colnames(x_tr)[which(coefs != 0)]
  }, error = function(e) {
    mspe_results["CR_Lasso"] <<- NA
    selected_vars[["CR_Lasso"]] <<- character(0)
  })
  
  # 6. RLARS (Proposed, K=1)
  tryCatch({
    fit <- srlars(x_tr, y_tr, n_models = 1, tolerance = 1e-8, 
                  x_preprocess = "ddc", y_preprocess = "wrap", cor_estimator = "wrap",
                  cv_preprocess = "global", cv_fit = "huber", cv_loss = "huber", compute_coef = TRUE)
    mspe_results["RLARS"] <- mean((predict(fit, newx = x_te) - y_te)^2)
    selected_vars[["RLARS"]] <- colnames(x_tr)[fit$active.sets[[1]]]
  }, error = function(e) {
    mspe_results["RLARS"] <<- NA
    selected_vars[["RLARS"]] <<- character(0)
  })
  
  # 7. FSCRE (Proposed, K=10)
  tryCatch({
    fit <- srlars(x_tr, y_tr, n_models = 10, tolerance = 1e-8, 
                  x_preprocess = "ddc", y_preprocess = "wrap", cor_estimator = "wrap",
                  cv_preprocess = "global", cv_fit = "huber", cv_loss = "huber", compute_coef = TRUE)
    mspe_results["FSCRE"] <- mean((predict(fit, newx = x_te) - y_te)^2)
    all_selected_idx <- unique(unlist(fit$active.sets))
    selected_vars[["FSCRE"]] <- colnames(x_tr)[all_selected_idx]
  }, error = function(e) {
    mspe_results["FSCRE"] <<- NA
    selected_vars[["FSCRE"]] <<- character(0)
  })
  
  names(mspe_results) <- paste0(prefix, "_", names(mspe_results))
  return(list(mspe = mspe_results, vars = selected_vars))
}


# ________________________________
# OUTER LOOP OVER TARGET PROTEINS
# ________________________________

proteins_to_run <- c("ER-alpha")
N_splits <- 50

for (target_protein in proteins_to_run) {
  
  cat(sprintf("\n\n*******************************************************\n"))
  cat(sprintf("STARTING ANALYSIS FOR TARGET PROTEIN: %s\n", target_protein))
  cat(sprintf("*******************************************************\n"))
  
  if (!(target_protein %in% colnames(prot_mat))) {
    warning("Target protein ", target_protein, " not found. Skipping.")
    next
  }
  
  y_full <- prot_mat[, target_protein]
  valid_y_idx <- which(!is.na(y_full))
  y <- y_full[valid_y_idx]
  X_full <- rna_mat[valid_y_idx, ]
  
  # Remove genes with near-zero variance globally just to drop empty columns
  gene_vars <- apply(X_full, 2, var)
  X_filtered <- X_full[, which(gene_vars > 1e-4)]
  
  mspe_list <- list()
  vars_list <- list() 
  
  for (i in 1:N_splits) {
    cat(sprintf("  Split %d / %d...\n", i, N_splits))
    
    # Train/Test Split
    train_idx <- sample(1:length(y), 50)
    
    x_tr_raw <- X_filtered[train_idx, ]
    y_tr_raw <- y[train_idx]
    x_te_raw <- X_filtered[-train_idx, ]
    y_te_raw <- y[-train_idx]
    
    # ________________________________________________
    # STRICT LEAKAGE-FREE FEATURE SCREENING & SCALING
    # ________________________________________________
    
    # 1. Correlation filter on Training Data ONLY
    gene_corrs <- abs(apply(x_tr_raw, 2, function(col) cor(col, y_tr_raw)))
    gene_corrs[is.na(gene_corrs)] <- 0
    top_genes <- order(gene_corrs, decreasing = TRUE)[1:500]
    
    x_train <- x_tr_raw[, top_genes]
    x_test  <- x_te_raw[, top_genes]
    
    # 2. Scale using Training Data parameters ONLY
    tr_means <- apply(x_train, 2, mean)
    tr_sds   <- apply(x_train, 2, sd)
    
    # Handle rare cases where SD is 0 after subsetting
    tr_sds[tr_sds == 0] <- 1 
    
    x_train <- scale(x_train, center = tr_means, scale = tr_sds)
    x_test  <- scale(x_test, center = tr_means, scale = tr_sds)
    
    y_mean <- mean(y_tr_raw)
    y_sd   <- sd(y_tr_raw)
    y_train <- as.numeric(scale(y_tr_raw, center = y_mean, scale = y_sd))
    y_test  <- as.numeric(scale(y_te_raw, center = y_mean, scale = y_sd))
    
    # A. Run on ORIGINAL Data
    res_orig <- runModels(x_train, y_train, x_test, y_test, prefix = "Orig")
    
    # B. Introduce TARGETED Artificial Contamination 
    x_train_cont <- x_train
    
    quick_en <- cv.glmnet(x_train, y_train, alpha = 3/4)
    clean_coefs <- as.numeric(coef(quick_en, s = "lambda.min"))[-1]
    
    n_top_genes <- min(30, sum(clean_coefs != 0))
    if (n_top_genes == 0) n_top_genes <- 30 
    top_var_indices <- order(abs(clean_coefs), decreasing = TRUE)[1:n_top_genes]
    
    n_to_corrupt <- round(0.15 * nrow(x_train))
    
    for (j in top_var_indices) {
        corrupt_rows <- sample(1:nrow(x_train), n_to_corrupt)
        x_train_cont[corrupt_rows, j] <- sample(c(10, -10), n_to_corrupt, replace = TRUE)
    }
    
    # Run on CONTAMINATED Data
    res_cont <- runModels(x_train_cont, y_train, x_test, y_test, prefix = "Contam")
    
    mspe_list[[i]] <- c(res_orig$mspe, res_cont$mspe)
    vars_list[[i]] <- list(Orig = res_orig$vars, Contam = res_cont$vars)
  }
  
  # ___________________________
  # Summarize and Save Results
  # ___________________________
  
  final_mspe_df <- do.call(rbind, mspe_list)
  
  cat(sprintf("\n--- FINAL RESULTS for %s (Average MSPE) ---\n", target_protein))
  print(round(colMeans(final_mspe_df, na.rm = TRUE), 4))
  
  safe_protein_name <- gsub("-", "_", target_protein)
  
  mspe_filename <- sprintf("results/Application_TCGA_%s_MSPE.rds", safe_protein_name)
  saveRDS(final_mspe_df, mspe_filename)
  
  vars_filename <- sprintf("results/Application_TCGA_%s_Genes.rds", safe_protein_name)
  saveRDS(vars_list, vars_filename)
  
  cat(sprintf("Results saved to %s and %s\n", mspe_filename, vars_filename))
}