#' @description generatePred() fits models and stores MSPE, recall, precision, and CPU time.
#' 
#' @param sim_data Output from generateData().
#' @param n_models Number of models for ensemble (K). 
#' @param tolerance Stopping tolerance for FSCRE.

# Required libraries
library(glmnet)
library(cellWise)
library(randomGLM)
library(parallel)
library(srlars) 
library(pense)
library(robustHD)
library(regcell) 

# Required source files
source("R/computeRCPR.R") 
source("SparseShootingS/sparseShootingS.R")

generatePred <- function(sim_data, n_models = 10, tolerance = 0, ...) {
  
  N <- length(sim_data$training_data$xtrain) 
  p.active <- sim_data$pactive 
  n <- sim_data$n
  p <- sim_data$p
  xtestdata <- sim_data$testing_data$xtest 
  ytestdata <- sim_data$testing_data$ytest
  true.beta <- sim_data$trueBeta
  
  # Define methods based on scenario
  if(sim_data$contamination.scenario != "casewise"){
    method_names <- c("ElasticNet", "DDC_EN", "DDC_RGLM", "Sparse_S", "CR_Lasso", "RLARS", "FSCRE")
  } else {
    method_names <- c("ElasticNet", "DDC_EN", "DDC_RGLM", "Sparse_S", "CR_Lasso", "RLARS", "FSCRE", "PENSE", "SparseLTS")
  }
  
  pred_output <- array(dim = c(length(method_names), 4, N))
  colnames(pred_output) <- c("MSPE","RC", "PR", "CPU")
  rownames(pred_output) <- method_names
  
  if(sim_data$contamination.scenario == "casewise") {
    cluster <- makeCluster(5)
  }
  
  for(i in 1:N) {
    cat("\n", "Iteration: ", i)
    
    xtrain <- sim_data$training_data$xtrain[[i]]
    ytrain <- sim_data$training_data$ytrain[[i]]
    
    # Ensure y is numeric vector
    ytrain <- as.numeric(ytrain)
    
    # ------------------------------------
    # 1. Baseline: Elastic Net (Raw Data)
    # ------------------------------------
    en_final <- tryCatch({
      cpu <- system.time(
        fit <- glmnet::cv.glmnet(x = xtrain, y = ytrain, alpha = 3/4)
      )["elapsed"]
      preds <- predict(fit, xtestdata, s = "lambda.min")
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["ElasticNet",, i] <- en_final
    
    # ------------------------------------------
    # 2. Shared Data Cleaning for DDC baselines
    # ------------------------------------------
    ddc_cpu <- 0
    x_imp <- xtrain
    y_imp <- ytrain
    
    tryCatch({
      ddc_cpu <- system.time({
        ddc_out <- cellWise::DDC(cbind(xtrain, ytrain), DDCpars = list(fastDDC = TRUE, silent = TRUE))
        x_imp <- ddc_out$Ximp[, 1:p]
        y_imp <- ddc_out$Ximp[, p+1]
      })["elapsed"]
    }, error = function(e) {
      warning("Shared DDC failed. Baselines will use raw data.")
    })

    # -------------------------------
    # 3. Baseline: DDC + Elastic Net
    # -------------------------------
    ddc_en_final <- tryCatch({
      cpu <- system.time(
        fit <- glmnet::cv.glmnet(x = x_imp, y = y_imp, alpha = 3/4)
      )["elapsed"]
      preds <- predict(fit, xtestdata, s = "lambda.min")
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu + ddc_cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["DDC_EN",, i] <- ddc_en_final
    
    # ------------------------------
    # 4. Baseline: DDC + Random GLM
    # ------------------------------
    ddc_rglm_final <- tryCatch({
      cpu <- system.time(
        # classify=FALSE is critical for regression
        fit <- randomGLM::randomGLM(x_imp, y_imp, classify = FALSE, nBags = 100, keepModels = TRUE, nThreads = 1, verbose = 0)
      )["elapsed"]
      preds <- predict(fit, newdata = xtestdata)
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      
      sel_counts <- colSums(fit$timesSelectedByForwardRegression)
      coefs <- rep(0, p)
      coefs[sel_counts > 0] <- 1
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu + ddc_cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["DDC_RGLM",, i] <- ddc_rglm_final

    # ---------------------------------------
    # 5. State-of-the-Art: Sparse Shooting S
    # ---------------------------------------
    sps_final <- tryCatch({
      cpu <- system.time(
        fit <- sparseshooting(x = xtrain, y = ytrain, wvalue = 3, nlambda = 50) 
      )["elapsed"]
      preds <- fit$coef[1] + xtestdata %*% fit$coef[-1]
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- fit$coef[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["Sparse_S",, i] <- sps_final
    
    # ------------------------------
    # 6. State-of-the-Art: CR-Lasso
    # ------------------------------
    cr_lasso_final <- tryCatch({
      cpu <- system.time(
        fit <- regcell::sregcell_std(ytrain, xtrain)
      )["elapsed"]
      
      preds <- fit$intercept_hat + xtestdata %*% fit$betahat
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- fit$betahat
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["CR_Lasso",, i] <- cr_lasso_final

    # -------------------------
    # 7. Proposed: RLARS (K=1)
    # -------------------------
    rlars_final <- tryCatch({
      cpu <- system.time(
        fit <- srlars::srlars(xtrain, ytrain, n_models = 1, tolerance = tolerance, robust = TRUE, compute_coef = TRUE)
      )["elapsed"]
      preds <- predict(fit, newx = xtestdata, dynamic = FALSE) 
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- coef(fit)[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["RLARS",, i] <- rlars_final

    # --------------------------
    # 8. Proposed: FSCRE (K=10)
    # --------------------------
    fscre_final <- tryCatch({
      cpu <- system.time(
        fit <- srlars::srlars(xtrain, ytrain, n_models = n_models, tolerance = tolerance, robust = TRUE, compute_coef = TRUE)
      )["elapsed"]
      preds <- predict(fit, newx = xtestdata, dynamic = FALSE) 
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- coef(fit)[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["FSCRE",, i] <- fscre_final

    # ----------------------
    # Casewise Only Methods
    # ----------------------
    if(sim_data$contamination.scenario == "casewise"){
      
      # PENSE
      pense_final <- tryCatch({
        cpu <- system.time(
          fit <- pense::adapense_cv(x = xtrain, y = ytrain, alpha = 0.75, cv_k = 5, cl = cluster)
        )["elapsed"]
        preds <- predict(fit, xtestdata)
        mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
        coefs <- coef(fit)[-1]
        metrics <- computeRCPR(coefs, sim_data$active_ind)
        c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
      }, error = function(e) c(NA, NA, NA, NA))
      pred_output["PENSE",, i] <- pense_final
      
      # Sparse LTS
      slts_final <- tryCatch({
        lambda_max <- robustHD::lambda0(xtrain, ytrain)
        lambda_grid = rev(exp(seq(log(1e-2*lambda_max), log(lambda_max), length = 20)))
        cpu <- system.time(
          fit <- robustHD::sparseLTS(x = xtrain, y = c(ytrain), lambda = lambda_grid, mode = "lambda", cluster = cluster)
        )["elapsed"]
        preds <- predict(fit, xtestdata)
        mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
        coefs <- coef(fit)[-1]
        metrics <- computeRCPR(coefs, sim_data$active_ind)
        c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
      }, error = function(e) c(NA, NA, NA, NA)) 
      pred_output["SparseLTS",, i] <- slts_final
    }
  }
  
  if(sim_data$contamination.scenario == "casewise") stopCluster(cluster)
  
  return(pred_output)
}

