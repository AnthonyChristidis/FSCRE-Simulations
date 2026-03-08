# ---------------------------------------------------------------------------
# Application 1: Pharmacogenomics (GDSC2) - Bortezomib Response Prediction
# ---------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# --- 1. Load Required Libraries ---
library(PharmacoGx)
library(glmnet)
library(cellWise)
library(randomGLM)
library(caret)
library(srlars) 
library(regcell)

# Load competitor scripts (adjust paths as needed)
source("SparseShootingS/sparseShootingS.R")

set.seed(0)

# --- 2. Data Downloading and Processing ---
cat("\n--- Preparing GDSC2 Dataset ---\n")

# Create data directory if it doesn't exist
if (!dir.exists("data")) dir.create("data")

# Use the exact PSet Name from availablePSets()
pset_name <- "GDSC_2020(v2-8.2)"
# Use the filename that PharmacoGx will automatically create
pset_file <- file.path("data", paste0(pset_name, ".rds"))

if(!file.exists(pset_file)) {
  cat(sprintf("Downloading %s (this may take a few minutes)...\n", pset_name))
  # PharmacoGx saves it automatically; we just assign the returned object
  GDSC2 <- PharmacoGx::downloadPSet(pset_name, saveDir = "data") 
} else {
  cat("Loading GDSC PSet from local file...\n")
  GDSC2 <- readRDS(pset_file)
}

# Extract RNA-seq expression data safely
cat("Extracting RNA-seq data...\n")
rna_obj <- molecularProfiles(GDSC2, "rna")

# Extract the data matrix and the cell line names based on the object type
if (inherits(rna_obj, "SummarizedExperiment")) {
  rna_data <- t(SummarizedExperiment::assay(rna_obj))
  # FIX: Use 'sampleid' based on the pheno output
  cell_line_names <- as.character(SummarizedExperiment::colData(rna_obj)$sampleid)
} else if (inherits(rna_obj, "ExpressionSet")) {
  rna_data <- t(Biobase::exprs(rna_obj))
  cell_line_names <- as.character(Biobase::pData(rna_obj)$sampleid)
} else {
  # If it returned a matrix directly
  rna_data <- t(rna_obj)
  pheno <- PharmacoGx::phenoInfo(GDSC2, "rna")
  cell_line_names <- as.character(pheno$sampleid)
}

# Handle potential replicates (keep first unique instance of each cell line)
unique_idx <- !duplicated(cell_line_names)
rna_data <- rna_data[unique_idx, ]
rownames(rna_data) <- cell_line_names[unique_idx]

# Extract Drug Sensitivity (AAC - Area Above Curve)
cat("Extracting Drug Sensitivity...\n")
drug_name <- "Bortezomib"
sens_data <- summarizeSensitivityProfiles(GDSC2, sensitivity.measure = "aac_recomputed")

if(!(drug_name %in% rownames(sens_data))) {
    stop("Drug not found. Available drugs include: ", paste(head(rownames(sens_data)), collapse=", "))
}

y_full <- sens_data[drug_name, ]

# Find overlapping cell lines (samples)
common_cells <- intersect(rownames(rna_data), names(y_full))

cat(sprintf("Found %d cell lines with both RNA and drug data.\n", length(common_cells)))

y <- y_full[common_cells]
X <- rna_data[common_cells, ]

# Remove NA responses
valid_idx <- which(!is.na(y))
y <- y[valid_idx]
X <- X[valid_idx, ]

# Filter to top 1000 most variable genes to create p >> n setting
cat("Filtering to top 1000 most variable genes...\n")
gene_vars <- apply(X, 2, var)
top_genes <- order(gene_vars, decreasing = TRUE)[1:1000]
X <- X[, top_genes]

# Standardize predictors (Important for penalized regression)
X <- scale(X)

cat(sprintf("Final Data Dimensions: n = %d, p = %d\n", nrow(X), ncol(X)))

# --- 3. Define Model Evaluation Helper ---
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
  
  # Prepend prefix to names (Original vs Contaminated)
  names(mspe_results) <- paste0(prefix, "_", names(mspe_results))
  return(mspe_results)
}

# --- 4. Simulation Loop ---
N_splits <- 50
results_list <- list()

cat("\n--- Starting Data Application Splits ---\n")

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
  
  # B. Introduce Artificial Contamination to Training Data (5% marginal shift)
  x_train_cont <- x_train
  n_cells <- nrow(x_train) * ncol(x_train)
  contam_idx <- sample(1:n_cells, round(0.05 * n_cells))
  x_train_cont[contam_idx] <- runif(length(contam_idx), min = 5, max = 15) # Strong shift
  
  # Run on CONTAMINATED Data
  mspe_cont <- run_models(x_train_cont, y_train, x_test, y_test, prefix = "Contam")
  
  # Store results
  results_list[[i]] <- c(mspe_orig, mspe_cont)
}

# --- 5. Summarize and Save Results ---
final_results_df <- do.call(rbind, results_list)

cat("\n--- FINAL RESULTS (Average MSPE) ---\n")
print(round(colMeans(final_results_df, na.rm = TRUE), 4))

if (!dir.exists("results")) dir.create("results")
saveRDS(final_results_df, "results/Application_GDSC_Results.rds")
cat("\nResults saved to results/Application1_GDSC_Results.rds\n")