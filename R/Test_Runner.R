# ---------------------------------------------
# FAST VERIFICATION SCRIPT for FSCRE Pipeline
# ---------------------------------------------

# This script runs a miniature version of the full simulation 
# to verify that all functions, loops, and data saving mechanisms 
# work without errors across all scenarios.

# Clear workspace
rm(list = ls())

# 1. Load Required Source Files
# (Ensure these are in the same directory or adjust paths)
source("R/computeRCPR.R")
source("R/generateData.R")
source("R/generatePred.R")
source("R/simFunc.R")
source("R/generateOutput.R")

# 2. Set "Toy" Parameters for Fast Testing
N_test <- 5              # Only 5 replications to test the inner loop
n_test <- 50             # Sample size
p_test <- 500            # Dimension
m_test <- 2000           # Test set
p_active_test <- 50      # Only test one sparsity level
snr_test <- 1            # Only test one SNR
rho_test <- 0.8
rho_inactive_test <- 0.2
group_size_test <- 5
n_models_test <- 5   # Ensemble size
sim_tolerance <- 1e-3

# Scenarios to test
scenarios_to_test <- c("casewise", 
                       "cellwise_marginal", 
                       "cellwise_correlation",
                       "mixture_marginal", 
                       "mixture_correlation")[5]

# Ensure output directory exists
if (!dir.exists("test_results")) {
  dir.create("test_results")
}

cat("=======================================================\n")
cat("STARTING FAST VERIFICATION RUN\n")
cat(sprintf("N=%d, n=%d, p=%d, p.active=%d\n", N_test, n_test, p_test, p_active_test[1]))
cat("=======================================================\n")

# 3. Main Loop
for(scenario_val in scenarios_to_test) {
  
  # Define test contamination proportions
  if(scenario_val == "casewise") {
    contam_prop <- c(0.1) # Just test 10%
  } else if (scenario_val %in% c("cellwise_marginal", "cellwise_correlation")) {
    contam_prop <- c(0.05) # Just test 5%
  } else if(scenario_val %in% c("mixture_marginal", "mixture_correlation")) {
    contam_prop <- c(0.1, 0.05) # Fixed mixture rates
  }
  
  cat(sprintf("\n>>> Testing Scenario: %s <<<\n", toupper(scenario_val)))
  
  # File name for test output
  filename <- paste0("test_results/test_scenario=", scenario_val, ".rds")
  
  # Run the pipeline
  # We use tryCatch so if one scenario fails, it doesn't stop the whole test script.
  test_result <- tryCatch({
    
    results <- generateOutput(N = N_test, 
                              n = n_test, 
                              m = m_test, 
                              p = p_test, 
                              rho = rho_test, 
                              rho.inactive = rho_inactive_test,
                              p.active = p_active_test, 
                              group.size = group_size_test, 
                              snr = snr_test, 
                              contamination.prop = contam_prop, 
                              contamination.scenario = scenario_val,
                              seed = 0, 
                              n_models = n_models_test,
                              tolerance = sim_tolerance)
    
    # Save the output
    saveRDS(results, file = filename)
    cat("  [SUCCESS] Scenario completed and saved to", filename, "\n")
    
    # Return true for success
    TRUE
    
  }, error = function(e) {
    cat("  [FAILED] Error in scenario:", e$message, "\n")
    return(FALSE)
  })
}

cat("\n=======================================================\n")
cat("VERIFICATION RUN COMPLETE.\n")
cat("Check 'test_results' folder for output files.\n")
cat("=======================================================\n")