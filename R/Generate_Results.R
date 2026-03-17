# ----------------------------------------
# Generate Results for FSCRE Simulations
# ----------------------------------------

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
rho <- 0.8

# Grid values
contamination_scenario <- c("casewise", "cellwise_marginal", "cellwise_correlation",
                            "mixture_marginal", "mixture_correlation")
snr_vec <- c(0.5, 1, 2)
p_active_vec <- c(50, 100, 200)

# Method Parameters
n_models <- 10        # K=10 as stated in Section 5.2
sim_tolerance <- 1e-5 # Low tolerance needed for this hard setting

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# -------------------------
# 2. Main Simulation Loop
# -------------------------

for(scenario_val in contamination_scenario) {
  
  # Determine the specific contamination proportions to loop over based on scenario
  if(scenario_val == "casewise") {
    props_to_test <- list(0, 0.1, 0.2) 
  } else if (scenario_val %in% c("cellwise_marginal", "cellwise_correlation")) {
    props_to_test <- list(0.05, 0.1) 
  } else if(scenario_val %in% c("mixture_marginal", "mixture_correlation")) {
    # For mixtures, it's a fixed pair of proportions, so just one item in the list
    props_to_test <- list(c(0.1, 0.05)) 
  }
    
  for(snr_val in snr_vec) {
    for(p_active_val in p_active_vec) {
      for(contam_val in props_to_test) {
        
        # Format the contamination string for the filename
        if (length(contam_val) > 1) {
          contam_str <- paste0(contam_val, collapse="_")
        } else {
          contam_str <- as.character(contam_val)
        }
        
        # Construct the highly specific filename
        filename <- paste0("results/res_scen=", scenario_val, 
                           "_snr=", snr_val, 
                           "_pAct=", p_active_val,
                           "_contam=", contam_str, ".rds")
        
        # Skip if already computed (great for restarting interrupted runs!)
        if (file.exists(filename)) {
          cat("\n Skipping: File already exists -", filename)
          next
        }
        
        # Generate results for this specific configuration
        results <- generateOutput(N = N, n = n, m = m, p = p, 
                                  rho = rho, rho.inactive = rho.inactive,
                                  p.active = p_active_val, group.size = group.size, 
                                  snr = snr_val, 
                                  contamination.prop = contam_val, 
                                  contamination.scenario = scenario_val,
                                  seed = 0, # Consider a dynamic seed if needed
                                  n_models = n_models,
                                  tolerance = sim_tolerance) 
        
        # Save JUST the single array object
        saveRDS(results, file = filename)
        cat("\n Results saved to:", filename, "\n")
        
      } # End contam loop
    } # End p.active loop
  } # End snr loop
} # End scenario loop

cat("\nAll simulations completed.\n")