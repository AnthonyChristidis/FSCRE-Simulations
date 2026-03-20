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

# Helper modified to return BOTH MSPE and selected variables
runModels <- function(x_tr, y_tr, x_te, y_te, prefix) {
  
  p <- ncol(x_tr)
  mspe_results <- c()
  selected_vars <- list()
  time_limit <- 300 # 5 minutes in seconds
  
  # 1. Elastic Net (EN)
  tryCatch({
    withTimeout({
      fit <- cv.glmnet(x = x_tr, y = y_tr, alpha = 3/4)
      mspe_results["ElasticNet"] <- mean((predict(fit, x_te, s = "lambda.min") - y_te)^2)
      
      # Extract non-zero coefficients (excluding intercept)
      coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
      selected_vars[["ElasticNet"]] <- colnames(x_tr)[which(coefs != 0)]
    }, timeout = time_limit, onTimeout = "error")
  }, error = function(e) {
    if (grepl("Timeout", e$message)) cat(" [EN Timed out] ")
    mspe_results["ElasticNet"] <<- NA
    selected_vars[["ElasticNet"]] <<- character(0)
  })
  
  # Shared DDC Imputation for Baselines
  x_imp <- x_tr
  y_imp <- y_tr
  tryCatch({
    withTimeout({
      ddc_out <- cellWise::DDC(cbind(x_tr, y_tr), DDCpars = list(fastDDC = TRUE, silent = TRUE))
      x_imp <<- ddc_out$Ximp[, 1:p]   
      y_imp <<- ddc_out$Ximp[, p+1]
    }, timeout = time_limit, onTimeout = "error")
  }, error = function(e) {
    if (grepl("Timeout", e$message)) cat(" [DDC Timed out] ")
  })
  
  # 2. DDC + EN
  tryCatch({
    withTimeout({
      fit <- cv.glmnet(x = x_imp, y = y_imp, alpha = 3/4)
      mspe_results["DDC_EN"] <- mean((predict(fit, x_te, s = "lambda.min") - y_te)^2)
      
      coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
      selected_vars[["DDC_EN"]] <- colnames(x_imp)[which(coefs != 0)]
    }, timeout = time_limit, onTimeout = "error")
  }, error = function(e) {
    if (grepl("Timeout", e$message)) cat(" [DDC_EN Timed out] ")
    mspe_results["DDC_EN"] <<- NA
    selected_vars[["DDC_EN"]] <<- character(0)
  })
  
  # 3. DDC + Random GLM
  tryCatch({
    withTimeout({
      fit <- randomGLM(x_imp, y_imp, classify = FALSE, nBags = 100, keepModels = TRUE, nThreads = 1, verbose = 0)
      mspe_results["DDC_RGLM"] <- mean((predict(fit, newdata = x_te) - y_te)^2)
      
      # Extract variables selected in > 0 bags
      sel_counts <- colSums(fit$timesSelectedByForwardRegression)
      selected_vars[["DDC_RGLM"]] <- colnames(x_imp)[which(sel_counts > 0)]
    }, timeout = time_limit, onTimeout = "error")
  }, error = function(e) {
    if (grepl("Timeout", e$message)) cat(" [DDC_RGLM Timed out] ")
    mspe_results["DDC_RGLM"] <<- NA
    selected_vars[["DDC_RGLM"]] <<- character(0)
  })
  
  # 4. Sparse Shooting S
  tryCatch({
    withTimeout({
      fit <- sparseshooting(x = x_tr, y = y_tr, wvalue = 3, nlambda = 50)
      preds <- fit$coef[1] + x_te %*% fit$coef[-1]
      mspe_results["Sparse_S"] <- mean((preds - y_te)^2)
      
      coefs <- fit$coef[-1]
      selected_vars[["Sparse_S"]] <- colnames(x_tr)[which(coefs != 0)]
    }, timeout = time_limit, onTimeout = "error")
  }, error = function(e) {
    if (grepl("Timeout", e$message)) cat(" [Sparse_S Timed out] ")
    mspe_results["Sparse_S"] <<- NA
    selected_vars[["Sparse_S"]] <<- character(0)
  })
  
  # 5. CR-Lasso
  tryCatch({
    withTimeout({
      fit <- regcell::sregcell_std(y_tr, x_tr)
      preds <- fit$intercept_hat + x_te %*% fit$betahat
      mspe_results["CR_Lasso"] <- mean((preds - y_te)^2)
      
      coefs <- fit$betahat
      selected_vars[["CR_Lasso"]] <- colnames(x_tr)[which(coefs != 0)]
    }, timeout = time_limit, onTimeout = "error")
  }, error = function(e) {
    if (grepl("Timeout", e$message)) cat(" [CR_Lasso Timed out] ")
    mspe_results["CR_Lasso"] <<- NA
    selected_vars[["CR_Lasso"]] <<- character(0)
  })
  
  # 6. RLARS (Proposed, K=1)
  tryCatch({
    withTimeout({
      fit <- srlars(x_tr, y_tr, n_models = 1, tolerance = 0.01, robust = TRUE, compute_coef = TRUE)
      mspe_results["RLARS"] <- mean((predict(fit, newx = x_te, dynamic = FALSE) - y_te)^2)
      
      # For K=1, active.sets is a list of length 1
      selected_vars[["RLARS"]] <- colnames(x_tr)[fit$active.sets[[1]]]
    }, timeout = time_limit, onTimeout = "error")
  }, error = function(e) {
    if (grepl("Timeout", e$message)) cat(" [RLARS Timed out] ")
    mspe_results["RLARS"] <<- NA
    selected_vars[["RLARS"]] <<- character(0)
  })
  
  # 7. FSCRE (Proposed, K=10)
  tryCatch({
    withTimeout({
      fit <- srlars(x_tr, y_tr, n_models = 10, tolerance = 0.01, robust = TRUE, compute_coef = TRUE)
      mspe_results["FSCRE"] <- mean((predict(fit, newx = x_te, dynamic = FALSE) - y_te)^2)
      
      # For K=10, active.sets is a list of K vectors. We want all unique genes selected across the ensemble.
      all_selected_idx <- unique(unlist(fit$active.sets))
      selected_vars[["FSCRE"]] <- colnames(x_tr)[all_selected_idx]
    }, timeout = time_limit, onTimeout = "error")
  }, error = function(e) {
    if (grepl("Timeout", e$message)) cat(" [FSCRE Timed out] ")
    mspe_results["FSCRE"] <<- NA
    selected_vars[["FSCRE"]] <<- character(0)
  })
  
  # Prepend prefix to names for MSPE
  names(mspe_results) <- paste0(prefix, "_", names(mspe_results))
  
  # We return a list containing both the MSPE vector and the list of selected variables
  return(list(mspe = mspe_results, vars = selected_vars))
}


# ________________________________
# OUTER LOOP OVER TARGET PROTEINS
# ________________________________

proteins_to_run <- c("EGFR", "ER-alpha")
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
  
  mspe_list <- list()
  
  # We will store the selected variables across all splits in a nested list
  # Structure: vars_list[[split_idx]][["Orig" or "Contam"]][["MethodName"]]
  vars_list <- list() 
  
  for (i in 1:N_splits) {
    cat(sprintf("  Split %d / %d...\n", i, N_splits))
    
    # Train/Test Split
    # Note: Using a fixed n_train = 50 to maintain the extreme p >> n setting as discussed
    train_idx <- sample(1:nrow(y), 50)
    
    x_train <- X[train_idx, ]
    y_train <- y[train_idx]
    x_test <- X[-train_idx, ]
    y_test <- y[-train_idx]
    
    # A. Run on ORIGINAL Data
    res_orig <- runModels(x_train, y_train, x_test, y_test, prefix = "Orig")
    
    # B. Introduce TARGETED Artificial Contamination 
    x_train_cont <- x_train
    
    # 1. Identify the "most important" variables from a quick, clean EN fit
    quick_en <- cv.glmnet(x_train, y_train, alpha = 3/4)
    clean_coefs <- as.numeric(coef(quick_en, s = "lambda.min"))[-1]
    
    # Get the indices of the top 30 most important genes
    n_top_genes <- min(30, sum(clean_coefs != 0))
    if (n_top_genes == 0) n_top_genes <- 30 
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
    res_cont <- runModels(x_train_cont, y_train, x_test, y_test, prefix = "Contam")
    
    # Store MSPE results
    mspe_list[[i]] <- c(res_orig$mspe, res_cont$mspe)
    
    # Store Variable Selection results
    vars_list[[i]] <- list(Orig = res_orig$vars, Contam = res_cont$vars)
  }
  
  # ___________________________
  # Summarize and Save Results
  # ___________________________
  
  final_mspe_df <- do.call(rbind, mspe_list)
  
  cat(sprintf("\n--- FINAL RESULTS for %s (Average MSPE) ---\n", target_protein))
  print(round(colMeans(final_mspe_df, na.rm = TRUE), 4))
  
  # Clean filename
  safe_protein_name <- gsub("-", "_", target_protein)
  
  # Save MSPE
  mspe_filename <- sprintf("results/Application_TCGA_%s_MSPE.rds", safe_protein_name)
  saveRDS(final_mspe_df, mspe_filename)
  
  # Save Variables
  vars_filename <- sprintf("results/Application_TCGA_%s_Genes.rds", safe_protein_name)
  saveRDS(vars_list, vars_filename)
  
  cat(sprintf("Results saved to %s and %s\n", mspe_filename, vars_filename))
}