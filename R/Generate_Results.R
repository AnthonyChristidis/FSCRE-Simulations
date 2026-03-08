#' @description Simulation Master Runner.

# Clear all memory
rm(list = ls())

# Required libraries
library(parallel)

# Required source file
source("R/generateOutput.R")

# ---------------------------------
# 1. Setting Simulation Parameters 
# ---------------------------------
N <- 50             # Replications
n <- 50             # Training sample size
p <- 500            # Total predictors
m <- 5000           # Test sample size
group.size <- 25    # Size of correlated blocks
rho.inactive <- 0.2 # Background correlation

# Grid values
contamination_scenario <- c("cellwise_marginal", "cellwise_correlation",
                            "mixture_marginal", "mixture_correlation",
                            "casewise")
snr <- c(0.5, 1, 2)
rho <- 0.8
p.active <- c(50, 100, 200)

# Method Parameters
n_models <- 10      # K=10 as stated in Section 5.2
sim_tolerance <- 1e-5 # Low tolerance needed for this hard setting

# Create results directory if it doesn't exist
if (!dir.exists("../results")) {
  dir.create("../results")
}

# -------------------------
# 2. Main Simulation Loop
# -------------------------
for(scenario_val in contamination_scenario) {
  
  # Define contamination proportions based on scenario
  if(scenario_val == "casewise") {
    contamination.prop <- c(0, 0.1, 0.2) 
  } else if (scenario_val %in% c("cellwise_marginal", "cellwise_correlation")) {
    contamination.prop <- c(0.05, 0.1) 
  } else if(scenario_val %in% c("mixture_marginal", "mixture_correlation")) {
    # e.g., 10% casewise, 5% cellwise
    contamination.prop <- c(0.1, 0.05) 
  }
    
  for(snr_val in snr) {
    for(rho_val in rho) {
      
      # Print simulation information
      cat("\n=======================================================")
      cat("\n Starting Batch:")
      cat("\n Contamination Scenario: ", scenario_val)
      cat("\n SNR: ", snr_val)
      cat("\n=======================================================\n")
      
      # File name with specifications of simulation (using .rds)
      filename <- paste0("../results/results_n=", n, "_p=", p, 
                         "_scenario=", scenario_val, 
                         "_snr=", snr_val, ".rds")
      
      # Generating results of simulation
      # Note: passing tolerance down the chain
      results <- generateOutput(N = N, n = n, m = m, p = p, 
                                rho = rho_val, rho.inactive = rho.inactive,
                                p.active = p.active, group.size = group.size, 
                                snr = snr_val, 
                                contamination.prop = contamination.prop, 
                                contamination.scenario = scenario_val,
                                seed = 0, 
                                n_models = n_models,
                                tolerance = sim_tolerance) # Pass tolerance
      
      # Saving JUST the simulation results object
      saveRDS(results, file = filename)
      cat("\n Results saved to:", filename, "\n")
    }
  }
}

# ----------------------------------------
# 3. Helper Function for Post-Processing
# ----------------------------------------
#' @description get_specific() extracts the output for MSPE, RC and PR for a specific setting.
get_specific <- function(result_obj,
                         target_p_active, 
                         target_contam_prop,
                         vec_p_active,
                         vec_contam_prop) {
  
  # Handle the case where contamination prop is a vector (mixture models)
  # In your current logic, mixture models return a flat list indexed by p.active
  if (is.list(target_contam_prop) || length(target_contam_prop) > 1) {
    p_idx <- which(vec_p_active == target_p_active)
    return(result_obj[[p_idx]])
  } else {
    c_idx <- which(vec_contam_prop == target_contam_prop)
    p_idx <- which(vec_p_active == target_p_active)
    return(result_obj[[c_idx]][[p_idx]])
  }
}