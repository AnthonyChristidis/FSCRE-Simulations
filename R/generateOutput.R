#' @description generateOutput() calls simFunc() for a single specific configuration.
#' 
#' @param N Number of training sets.
#' @param n Sample size for training set.
#' @param m Sample size for test set.
#' @param p Total number of parameters.
#' @param rho Correlation within a block of active parameters.
#' @param rho.inactive Correlation between blocks of active parameters.
#' @param p.active A single active parameter count.
#' @param group.size Size of one block of active parameters.
#' @param snr Signal to noise ratio.
#' @param contamination.prop A single contamination proportion (or a specific vector of 2 for mixture).
#' @param contamination.scenario Scenario string.
#' @param seed Random seed.
#' @param n_models Number of models for ensemble.

# Required source file
source("R/simFunc.R")

generateOutput <- function (N, 
                            n, 
                            m, 
                            p, 
                            rho, 
                            rho.inactive = 0.2,
                            p.active, 
                            group.size,
                            snr, 
                            contamination.prop, 
                            contamination.scenario,
                            seed = 0,
                            n_models,
                            ...){
  # Setting seed
  set.seed(seed)
  
  # Print current settings for the log
  cat("\n============================================\n")
  cat("Scenario:", contamination.scenario, "\n")
  if (length(contamination.prop) == 1) {
    cat("Contamination prop:", contamination.prop, "\n")
  } else {
    cat("Contamination prop:", paste(contamination.prop, collapse=", "), "\n")
  }
  cat("Active predictors:", p.active, "\n")
  cat("SNR:", snr, "\n")
  cat("============================================\n")
  
  # Compute simFunc() for this single specific setting
  output <- simFunc(N = N, n = n, m = m, p = p, 
                    rho = rho, rho.inactive = rho.inactive, 
                    p.active = p.active, 
                    group.size = group.size, snr = snr, 
                    contamination.prop = contamination.prop, 
                    contamination.scenario = contamination.scenario,
                    n_models = n_models, 
                    ...)
  
  return(output)
}