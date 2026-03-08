#' @description simfunc() calls generateData() and generatePred().
#' 
#' @param N Number of training sets (replications).
#' @param n Sample size for training set.
#' @param m Sample size for test set.
#' @param p Total number of predictors.
#' @param rho Correlation within a block of active predictors.
#' @param rho.inactive Correlation between blocks of active predictors.
#' @param p.active Number of true active predictors.
#' @param group.size Size of one block of active predictors.
#' @param snr Signal-to-noise ratio.
#' @param contamination.prop Contamination proportion.
#' @param contamination.scenario String specifying the contamination type.
#' @param n_models Number of models for the FSCRE ensemble (K).
#' @param ... Additional arguments passed to generatePred (e.g., tolerance).

# Required source files 
source("R/generateData.R")
source("R/generatePred.R")

simfunc <- function(N, 
                    n, 
                    m, 
                    p, 
                    rho, 
                    rho.inactive,
                    p.active, 
                    group.size,
                    snr, 
                    contamination.prop,
                    contamination.scenario,
                    n_models,
                    ...){ 
  
  # 1. Generate the datasets
  sim_data <- generateData(N = N, n = n, m = m, p = p, 
                           rho = rho, rho.inactive = rho.inactive, 
                           p.active = p.active, group.size = group.size, 
                           snr = snr, 
                           contamination.prop = contamination.prop, 
                           contamination.scenario = contamination.scenario)
                           
  # 2. Run the methods and collect metrics
  output <- generatePred(sim_data = sim_data, n_models = n_models, ...)
  
  return(output)
}