# -------------------------------------------
# Computational Scalability Study
# -------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries
library(mvnfast)
library(glmnet)
library(cellWise)
library(srlars)

# -------------------------------------------

# Set seed for reproducibility
set.seed(0)

# ______________________
# Simulation Parameters
# ______________________

# Grids for n and p
n_grid <- c(50, 100, 200, 500, 1000)
p_grid <- c(50, 100, 200, 500, 1000, 2000, 5000, 10000)

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
sim_tolerance <- 1e-5

# Data structure to hold results
timing_results <- data.frame(
  n = integer(),
  p = integer(),
  Method = character(),
  Rep = integer(),
  Time = numeric()
)

# Ensure output directory exists
if (!dir.exists("results")) dir.create("results")

# _________________
# Main Timing Loop
# _________________

cat("\n--- Starting Computational Scalability Study ---\n")

for (n in n_grid) {
  for (p in p_grid) {
    
    # We enforce p.active to be min(50, p) so we don't ask for more active vars than exist
    p.active <- min(50, p) 
    
    cat(sprintf("\nRunning Setting: n = %d, p = %d (p.active = %d)\n", n, p, p.active))
    
    for (rep in 1:N_reps) {
      
      # _________________________________________
      # A. Data Generation (Mixture Correlation)
      # _________________________________________
      
      # Block Correlation
      sigma.mat <- matrix(0, nrow = p, ncol = p)
      sigma.mat[1:p.active, 1:p.active] <- rho.inactive
      if (p.active >= group.size) {
        for(group in 0:(p.active/group.size - 1)) {
          idx <- (group*group.size+1):(group*group.size+group.size)
          sigma.mat[idx, idx] <- rho
        }
      }
      diag(sigma.mat) <- 1
      
      # True Beta
      trueBeta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
      
      # Noise
      sigma_val <- as.numeric(sqrt(t(trueBeta) %*% sigma.mat %*% trueBeta)/sqrt(snr))
      
      # Clean Data
      x_train <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
      y_train <- as.numeric(x_train %*% trueBeta + rnorm(n, 0, sigma_val))
      
      # Contamination Setup
      beta_cont <- trueBeta
      beta_cont[trueBeta!=0] <- beta_cont[trueBeta!=0]*(1 + k_slo)
      beta_cont[trueBeta==0] <- k_slo*max(abs(trueBeta))
      if (max(abs(trueBeta)) == 0) beta_cont[trueBeta==0] <- k_slo # edge case safety
      
      n.casewise <- floor(n * contam_prop_casewise)
      
      # 1. Casewise Contamination
      if (n.casewise > 0) {
        contamination_indices <- 1:n.casewise
        for(cont_id in contamination_indices){
          a <- runif(p, min = -1, max = 1)
          a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
          x_train[cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
            k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
          y_train[cont_id] <- t(x_train[cont_id,]) %*% beta_cont
        }
      }
      
      # 2. Cellwise Correlation Contamination
      if (n - n.casewise > 0) {
        n_clean_rows <- n - n.casewise
        contamination_indices <- sample(1:(n_clean_rows * p), round(n_clean_rows * p * contam_prop_cellwise))
        subset_idx <- (n.casewise + 1):n
        
        sub_matrix <- x_train[subset_idx, , drop=FALSE]
        sub_matrix[contamination_indices] <- NA
        
        for(row_id in 1:nrow(sub_matrix)){
          cells_id <- which(is.na(sub_matrix[row_id,]))
          if(length(cells_id) > 0) {
            if (length(cells_id) > 1) {
              mu_cells <- rep(0, length(cells_id))
              sigma_cells <- sigma.mat[cells_id, cells_id, drop=FALSE]
              eigen_vec <- eigen(sigma_cells)$vectors[, length(cells_id)]
              sub_matrix[row_id, cells_id] <- gamma * sqrt(length(cells_id)) * t(eigen_vec) /
                sqrt(mahalanobis(t(eigen_vec), mu_cells, sigma_cells))
            } else {
              sub_matrix[row_id, cells_id] <- gamma * 3
            }
          }
        }
        x_train[subset_idx, ] <- sub_matrix
      }
      
      # _______________________
      # B. Timing FSCRE (K=10)
      # _______________________

      time_fscre <- system.time({
        # We suppress warnings to keep the console clean during the long loop
        suppressWarnings({
          fit_fscre <- srlars(x_train, y_train, 
                              n_models = n_models, 
                              tolerance = sim_tolerance, 
                              robust = TRUE, 
                              compute_coef = TRUE)
        })
      })["elapsed"]
      
      timing_results <- rbind(timing_results, data.frame(
        n = n, p = p, Method = "FSCRE", Rep = rep, Time = unname(time_fscre)
      ))
      
      # _____________________________
      # C. Timing DDC + Elastic Net
      # _____________________________
      
      time_ddc_en <- system.time({
        suppressWarnings({
          # Shared DDC Imputation
          ddc_out <- cellWise::DDC(cbind(x_train, y_train), DDCpars = list(fastDDC = TRUE, silent = TRUE))
          x_imp <- ddc_out$Ximp[, 1:p]
          y_imp <- ddc_out$Ximp[, p+1]
          
          # Elastic Net
          fit_en <- cv.glmnet(x = x_imp, y = y_imp, alpha = 0.5)
        })
      })["elapsed"]
      
      timing_results <- rbind(timing_results, data.frame(
        n = n, p = p, Method = "DDC_EN", Rep = rep, Time = unname(time_ddc_en)
      ))
      
    } # End Reps
    
    # Save intermediate results so we don't lose data if it crashes on p=10000
    saveRDS(timing_results, file = "results/Timing_Results_Intermediate.rds")
    
  } # End p_grid
} # End n_grid

# ___________
# Final Save
# ___________

cat("\n--- Timing Study Complete ---\n")
saveRDS(timing_results, file = "results/Timing_Results_Final.rds")
cat("Results saved to results/Timing_Results_Final.rds\n")