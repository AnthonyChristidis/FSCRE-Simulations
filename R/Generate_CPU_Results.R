# -------------------------------------------
# Computational Scalability Study
# -------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(mvnfast)
library(glmnet)
library(cellWise)
library(randomGLM)
library(srlars)

# -------------------------------------------

# Set seed for reproducibility
set.seed(0)

# ______________________
# Simulation Parameters
# ______________________

N_reps <- 50 # Number of replications per setting

# Fixed parameters for the timing study
snr <- 1
rho <- 0.8
rho.inactive <- 0.2
group.size <- 25
contamination.scenario <- "mixture_correlation"
contam_prop_casewise <- 0.10
contam_prop_cellwise <- 0.05
gamma <- 3
k_lev <- 2
k_slo <- 100

# FSCRE parameters
n_models <- 10
sim_tolerance <- 1e-8 # Strict, non-negative tolerance

# ___________________________
# Define Full Grid Settings
# ___________________________

n_grid <- c(50, 100, 200, 500)
p_grid <- c(50, 100, 500, 1000, 2000, 5000)

# Create a data frame of ALL combinations (4 * 6 = 24 unique settings)
settings <- expand.grid(n = n_grid, p = p_grid)

# Ensure output directory exists
if (!dir.exists("results")) {
  dir.create("results")
}

# _________________
# Main Timing Loop
# _________________

cat("\n--- Starting Computational Scalability Study ---\n")

for (i in 1:nrow(settings)) {
  
  n <- settings$n[i]
  p <- settings$p[i]
  p.active <- min(50, p) 
  
  filename <- paste0("results/timing_n=", n, "_p=", p, ".rds")
  
  if (file.exists(filename)) {
    cat(sprintf("\nSkipping Setting %d/%d: n=%d, p=%d\n", i, nrow(settings), n, p))
    next 
  }
  
  cat(sprintf("\nRunning Setting %d/%d: n = %d, p = %d (p.active = %d)\n", i, nrow(settings), n, p, p.active))
  
  # These do not change across replications, so compute them ONCE per grid point.
  
  # 1. Block Correlation Matrix
  sigma.mat <- matrix(rho.inactive, nrow = p, ncol = p)
  if (p.active >= group.size) {
    for(group in 0:(p.active/group.size - 1)) {
      idx <- (group*group.size+1):(group*group.size+group.size)
      sigma.mat[idx, idx] <- rho
    }
  }
  diag(sigma.mat) <- 1
  
  # 2. Inverse of Sigma (for Casewise leverage scaling)
  # solve() on p=5000 is slow. We do it ONCE here.
  inv_sigma <- solve(sigma.mat)
  
  # 3. Smallest Eigenvector of Global Sigma (for Cellwise Correlation)
  # Instead of computing sub-eigenvectors per row (O(n * p^3)), we compute the global 
  # smallest eigenvector once. Perturbing cells along this global direction of least 
  # variance is just as adversarial and orders of magnitude faster.
  eig_res <- eigen(sigma.mat, symmetric = TRUE)
  global_e_min <- eig_res$vectors[, p] # Last column is smallest eigenvalue
  
  # 4. True Beta
  trueBeta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
  sigma_val <- as.numeric(sqrt(t(trueBeta) %*% sigma.mat %*% trueBeta)/sqrt(snr))
  
  # 5. Distorted Beta (for Casewise)
  beta_cont <- trueBeta
  beta_cont[trueBeta!=0] <- beta_cont[trueBeta!=0]*(1 + k_slo)
  beta_cont[trueBeta==0] <- k_slo*max(abs(trueBeta))
  if (max(abs(trueBeta)) == 0) beta_cont[trueBeta==0] <- k_slo 
  
  n.casewise <- floor(n * contam_prop_casewise)
  n_clean_rows <- n - n.casewise
  
  setting_results_list <- list()
  counter <- 1
  
  for (rep in 1:N_reps) {

    cat("Rep:", rep, "\n")
    
    # _________________________________________
    # A. Data Generation (Optimized)
    # _________________________________________
    
    # Clean Data
    x_train <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
    colnames(x_train) <- paste0("V", 1:p) # Ensure colnames for DDC compatibility
    y_train <- as.numeric(x_train %*% trueBeta + rnorm(n, 0, sigma_val))
    
    # 1. Casewise Contamination (Vectorized)
    if (n.casewise > 0) {
      cont_idx <- 1:n.casewise
      
      # Generate all random 'a' vectors at once: matrix of size (n.casewise x p)
      A_mat <- matrix(runif(n.casewise * p, min = -1, max = 1), nrow = n.casewise, ncol = p)
      # Center rows
      A_mat <- t(apply(A_mat, 1, function(x) x - mean(x)))
      
      # Compute scale factor for each row: a^T * inv(Sigma) * a
      scale_factors <- sqrt(rowSums((A_mat %*% inv_sigma) * A_mat))
      
      # Generate multivariate normal noise for all contaminated rows at once
      noise_mat <- mvnfast::rmvn(n.casewise, mu = rep(0, p), sigma = 0.1^2*diag(p))
      
      # Apply contamination
      x_train[cont_idx, ] <- noise_mat + k_lev * (A_mat / scale_factors)
      y_train[cont_idx] <- as.numeric(x_train[cont_idx, ] %*% beta_cont)
    }
    
    # 2. Cellwise Correlation Contamination (Optimized)
    if (n_clean_rows > 0) {
      subset_idx <- (n.casewise + 1):n
      n_cells_to_contaminate <- round(n_clean_rows * p * contam_prop_cellwise)
      
      # 1. Randomly select cell indices to contaminate
      contam_rows <- sample(subset_idx, n_cells_to_contaminate, replace = TRUE)
      contam_cols <- sample(1:p, n_cells_to_contaminate, replace = TRUE)
      
      # 2. Apply adversarial shift based on global smallest eigenvector
      shifts <- gamma * 3 * global_e_min[contam_cols]
      
      # Apply the shifts vectorized
      x_train[cbind(contam_rows, contam_cols)] <- shifts
    }

    # _______________________
    # B. Timing FSCRE (K=10)
    # _______________________

    time_fscre <- system.time({
      suppressWarnings({
        fit_fscre <- srlars(x_train, y_train, 
                            n_models = n_models, 
                            tolerance = sim_tolerance, 
                            x_preprocess = "ddc",
                            y_preprocess = "wrap",
                            cor_estimator = "wrap",
                            cv_preprocess = "global",
                            cv_fit = "huber",
                            cv_loss = "huber",
                            compute_coef = TRUE)
      })
    })["elapsed"]
    
    setting_results_list[[counter]] <- data.frame(n = n, p = p, Method = "FSCRE", Rep = rep, Time = unname(time_fscre))
    counter <- counter + 1
    
    # ______________________________________
    # C. Timing Baseline Data Cleaning (DDC)
    # ______________________________________
    
    time_ddc <- system.time({
      suppressWarnings({
        # DDC Imputation on X ONLY
        ddc_out <- cellWise::DDC(x_train, DDCpars = list(fastDDC = TRUE, silent = TRUE))
        x_imp <- ddc_out$Ximp
        
        # Wrap response univariately
        y_imp <- as.numeric(cellWise::wrap(as.matrix(y_train))$Xw)
      })
    })["elapsed"]
    
    # _____________________________
    # D. Timing DDC + Elastic Net
    # _____________________________
    
    time_en <- system.time({
      suppressWarnings({
        fit_en <- glmnet::cv.glmnet(x = x_imp, y = y_imp, alpha = 0.5)
      })
    })["elapsed"]
    
    # Total pipeline time = DDC preprocessing + EN fitting
    total_time_ddc_en <- time_ddc + time_en
    setting_results_list[[counter]] <- data.frame(n = n, p = p, Method = "DDC_EN", Rep = rep, Time = unname(total_time_ddc_en))
    counter <- counter + 1
    
    # ___________________________
    # E. Timing DDC + Random GLM
    # ___________________________
    
    time_rglm <- system.time({
      suppressWarnings({
        fit_rglm <- randomGLM::randomGLM(x_imp, y_imp, 
                                         classify = FALSE, 
                                         nBags = 100, 
                                         keepModels = TRUE, 
                                         nThreads = 1, 
                                         verbose = 0)
      })
    })["elapsed"]
    
    # Total pipeline time = DDC preprocessing + RGLM fitting
    total_time_ddc_rglm <- time_ddc + time_rglm
    setting_results_list[[counter]] <- data.frame(n = n, p = p, Method = "DDC_RGLM", Rep = rep, Time = unname(total_time_ddc_rglm))
    counter <- counter + 1
    
  } # End Reps
  
  # _________________________________________
  # F. Save Results for this (n, p) setting
  # _________________________________________
  
  setting_results_df <- do.call(rbind, setting_results_list)
  saveRDS(setting_results_df, file = filename)
  cat(sprintf("  Saved %s\n", filename))
  
} # End settings loop

cat("\n--- Timing Study Complete ---\n")