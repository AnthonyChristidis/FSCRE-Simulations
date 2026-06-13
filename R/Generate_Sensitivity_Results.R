# -------------------------------------------------------------
# Sensitivity Analysis for FSCRE Number of Models (K)
# -------------------------------------------------------------

# Clear all memory
rm(list = ls())

# Required libraries
library(srlars)
library(ggplot2)
library(dplyr)
library(tidyr)

# Required source files
source("R/generateData.R")
source("R/computeRCPR.R") 

# -------------------------------------------------------------

# _________________________________________________
# 1. Custom Prediction Function for K-Sensitivity
# _________________________________________________

generatePred_K_Sensitivity <- function(sim_data, K_vec, tolerance = 1e-8) {
  
  N <- length(sim_data$training_data$xtrain) 
  p <- sim_data$p
  xtestdata <- sim_data$testing_data$xtest 
  ytestdata <- sim_data$testing_data$ytest
  
  # Ensure column names for DDCpredict compatibility
  colnames(xtestdata) <- paste0("V", 1:p)
  
  # Setup output array: Rows = K values, Cols = Metrics, Slices = Replications
  pred_output <- array(dim = c(length(K_vec), 4, N))
  colnames(pred_output) <- c("MSPE", "RC", "PR", "CPU")
  rownames(pred_output) <- paste0("K_", K_vec)
  
  for(i in 1:N) {
    cat("\n Iteration:", i, "of", N, "| Evaluating K =")
    
    xtrain <- sim_data$training_data$xtrain[[i]]
    colnames(xtrain) <- paste0("V", 1:p)
    ytrain <- as.numeric(sim_data$training_data$ytrain[[i]])
    
    for(k_idx in seq_along(K_vec)) {
      K <- K_vec[k_idx]
      cat(" ", K)
      
      # _______________________________
      # Evaluate FSCRE for current K
      # _______________________________

      fscre_final <- tryCatch({
        cpu <- system.time({
          fit <- srlars::srlars(xtrain, ytrain, 
                                n_models = K, 
                                tolerance = tolerance, 
                                x_preprocess = "ddc",
                                y_preprocess = "wrap",
                                cor_estimator = "wrap",
                                cv_preprocess = "global",
                                cv_fit = "huber",
                                cv_loss = "huber",
                                compute_coef = TRUE)
          # Predict automatically handles DDCpredict on the test set
          preds <- predict(fit, newx = xtestdata) 
        })["elapsed"]
        
        mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
        coefs <- coef(fit)[-1]
        metrics <- computeRCPR(coefs, sim_data$active_ind)
        
        c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
      }, error = function(e) c(NA, NA, NA, NA))
      
      pred_output[k_idx, , i] <- fscre_final
    }
  }
  cat("\n")
  return(pred_output)
}

# __________________________________
# 2. Setting Simulation Parameters 
# __________________________________

N <- 50             # Replications
n <- 50             # Training sample size
p <- 500            # Total predictors
m <- 5000           # Test sample size
group.size <- 25    # Size of correlated blocks
rho.inactive <- 0.2 # Background correlation
rho <- 0.8

# Fixed parameters for this specific sensitivity study
scenario_val <- "mixture_correlation"
snr_val <- 1.0            # Moderate signal
contam_val <- c(0.1, 0.05)
sim_tolerance <- 1e-8     # Strict, non-negative tolerance

# Loop grids
p_active_vec <- c(50, 100, 200)             # High, Moderate, and Low sparsity
K_vec <- as.numeric(1:20)                   # Evaluate K from 1 to 20

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# ________________________
# 3. Main Simulation Loop
# ________________________

# Formatting the contamination string for the filename
contam_str <- paste0(contam_val, collapse="_")

for(p_active_val in p_active_vec) {
  
  # Construct the highly specific filename
  filename <- paste0("results/sensitivity_K_scen=", scenario_val, 
                     "_snr=", snr_val, 
                     "_pAct=", p_active_val,
                     "_contam=", contam_str, ".rds")
  
  # Skip if already computed 
  if (file.exists(filename)) {
    cat("\n Skipping: File already exists -", filename)
    next
  }
  
  cat("\n============================================\n")
  cat("Sensitivity Analysis for K (FSCRE)\n")
  cat("Scenario:", scenario_val, "\n")
  cat("Active predictors:", p_active_val, "\n")
  cat("SNR:", snr_val, "\n")
  cat("============================================\n")
  
  # 1. Generate the datasets
  set.seed(0) # Keep seed consistent with main simulations
  sim_data <- generateData(N = N, n = n, m = m, p = p, 
                           rho = rho, rho.inactive = rho.inactive, 
                           p.active = p_active_val, group.size = group.size, 
                           snr = snr_val, 
                           contamination.prop = contam_val, 
                           contamination.scenario = scenario_val)
  
  # 2. Run the custom K-sensitivity method
  results <- generatePred_K_Sensitivity(sim_data = sim_data, 
                                        K_vec = K_vec, 
                                        tolerance = sim_tolerance)
  
  # 3. Save JUST the single array object
  saveRDS(results, file = filename)
  cat("Results saved to:", filename, "\n")
  
} # End p.active loop

cat("\nAll sensitivity simulations completed.\n")

# _____________________________________________________
# 4. Process Results and Generate Presentation Figure
# _____________________________________________________

# Parameters used in the simulation
scenario_val <- "mixture_correlation"
snr_val <- 1.0
contam_str <- "0.1_0.05"
p_active_vec <- c(50, 100, 200)
K_vec <- 1:20

# Initialize an empty data frame to store aggregated results
summary_data <- data.frame()

# Loop through the saved files and extract median performance
for (p_active_val in p_active_vec) {
  
  filename <- paste0("results/sensitivity_K_scen=", scenario_val, 
                     "_snr=", snr_val, 
                     "_pAct=", p_active_val,
                     "_contam=", contam_str, ".rds")
  
  if (file.exists(filename)) {
    res_array <- readRDS(filename)
    res_median <- apply(res_array, c(1, 2), median, na.rm = TRUE)
    
    df_temp <- as.data.frame(res_median)
    df_temp$K <- K_vec
    df_temp$p_active <- paste0("Active Predictors: ", p_active_val, " (", (p_active_val/500)*100, "%)")
    
    summary_data <- bind_rows(summary_data, df_temp)
  } else {
    warning(paste("File not found:", filename))
  }
}

# _____________________________
# Data Formatting for ggplot2
# _____________________________

plot_data <- summary_data |>
  select(K, p_active, MSPE, RC, PR) |>
  pivot_longer(cols = c(MSPE, RC, PR), names_to = "Metric", values_to = "Value")

plot_data$Metric <- factor(plot_data$Metric, 
                           levels = c("MSPE", "RC", "PR"), 
                           labels = c("Prediction Error (MSPE)", "Recall", "Precision"))

plot_data$p_active <- factor(plot_data$p_active, 
                             levels = c("Active Predictors: 50 (10%)", 
                                        "Active Predictors: 100 (20%)", 
                                        "Active Predictors: 200 (40%)"))

# Extract the RLARS Baseline (K=1) to plot as a horizontal line
rlars_data <- plot_data |> filter(K == 1)

# 
# --- PUBLICATION-READY THEME & COLORS ---
# 

pub_theme <- theme_bw() +
  theme(
    text = element_text(size = 14, family = "sans"),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12, color = "black"),
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 13)
  )

fscre_color <- "#0072B2"
rlars_color <- "#D55E00"

# _____________________
# Generate the Figure
# _____________________

if (!dir.exists("figures")) dir.create("figures")

p <- ggplot(plot_data, aes(x = K, y = Value)) +
  # Add the RLARS (K=1) horizontal baseline
  geom_hline(data = rlars_data, aes(yintercept = Value, linetype = "RLARS (Single Model)"), 
             color = rlars_color, linewidth = 1) +
  # Add the FSCRE curve
  geom_line(aes(color = "FSCRE (Ensemble)"), linewidth = 1.2) +
  geom_point(aes(color = "FSCRE (Ensemble)"), size = 2.5) +
  facet_grid(Metric ~ p_active, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_color_manual(values = c("FSCRE (Ensemble)" = fscre_color)) +
  scale_linetype_manual(values = c("RLARS (Single Model)" = "dashed")) +
  pub_theme +
  labs(
    x = expression(bold("Number of Models (") * italic(K) * bold(")")),
    y = "Median Performance Value"
  )

print(p)

ggsave(filename = "figures/FSCRE_K_Sensitivity.pdf", plot = p, 
       width = 10, height = 7, units = "in", dpi = 300)

cat("\nFigure saved to: figures/FSCRE_K_Sensitivity.pdf\n")