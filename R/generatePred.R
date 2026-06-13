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

generatePred <- function(sim_data, n_models = 10, tolerance = 1e-8, ...) {
  
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
    
    # Add column names to ensure DDCpredict works properly across train and test
    colnames(xtrain) <- paste0("V", 1:p)
    colnames(xtestdata) <- paste0("V", 1:p)
    
    # ____________________________________
    # 1. Baseline: Elastic Net (Raw Data)
    # ____________________________________

    en_final <- tryCatch({
      cpu <- system.time({
        fit <- glmnet::cv.glmnet(x = xtrain, y = ytrain, alpha = 3/4)
        preds <- predict(fit, xtestdata, s = "lambda.min")
      })["elapsed"]
      
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["ElasticNet",, i] <- en_final
    
    # _______________________________________________________________
    # 2. Shared Data Cleaning for DDC baselines (X-Only, No Leakage)
    # _______________________________________________________________

    ddc_cpu <- 0
    x_imp <- xtrain
    y_imp <- ytrain
    xtest_imp <- xtestdata
    
    tryCatch({
      ddc_cpu <- system.time({
        # Clean predictors ONLY
        ddc_out <- cellWise::DDC(xtrain, DDCpars = list(fastDDC = TRUE, silent = TRUE))
        x_imp <- ddc_out$Ximp
        
        # Apply the trained DDC to the test set predictors
        ddc_test <- cellWise::DDCpredict(xtestdata, ddc_out)
        xtest_imp <- ddc_test$Ximp
        
        # Robust univariate wrap for the response
        y_imp <- as.numeric(cellWise::wrap(as.matrix(ytrain))$Xw)
      })["elapsed"]
    }, error = function(e) {
      warning("Shared DDC failed. Baselines will use raw data.")
    })

    # ________________________________
    # 3. Baseline: DDC + Elastic Net
    # ________________________________

    ddc_en_final <- tryCatch({
      cpu <- system.time({
        fit <- glmnet::cv.glmnet(x = x_imp, y = y_imp, alpha = 3/4)
        # Predict on the DDC-cleaned test data
        preds <- predict(fit, xtest_imp, s = "lambda.min")
      })["elapsed"]
      
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- as.numeric(coef(fit, s = "lambda.min"))[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu + ddc_cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["DDC_EN",, i] <- ddc_en_final
    
    # _______________________________
    # 4. Baseline: DDC + Random GLM
    # _______________________________

    ddc_rglm_final <- tryCatch({
      cpu <- system.time({
        fit <- randomGLM::randomGLM(x_imp, y_imp, classify = FALSE, nBags = 100, keepModels = TRUE, nThreads = 1, verbose = 0)
        # Predict on the DDC-cleaned test data
        preds <- predict(fit, newdata = xtest_imp)
      })["elapsed"]
      
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      sel_counts <- colSums(fit$timesSelectedByForwardRegression)
      coefs <- rep(0, p)
      coefs[sel_counts > 0] <- 1
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu + ddc_cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["DDC_RGLM",, i] <- ddc_rglm_final

    # _______________________________________
    # 5. State-of-the-Art: Sparse Shooting S
    # _______________________________________

    sps_final <- tryCatch({
      cpu <- system.time({
        fit <- sparseshooting(x = xtrain, y = ytrain, wvalue = 3, nlambda = 50) 
        preds <- fit$coef[1] + xtestdata %*% fit$coef[-1]
      })["elapsed"]
      
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- fit$coef[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["Sparse_S",, i] <- sps_final
    
    # ______________________________
    # 6. State-of-the-Art: CR-Lasso
    # ______________________________

    cr_lasso_final <- tryCatch({
      cpu <- system.time({
        fit <- regcell::sregcell_std(ytrain, xtrain)
        preds <- fit$intercept_hat + xtestdata %*% fit$betahat
      })["elapsed"]
      
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- fit$betahat
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["CR_Lasso",, i] <- cr_lasso_final

    # _________________________
    # 7. Proposed: RLARS (K=1)
    # _________________________

    rlars_final <- tryCatch({
      cpu <- system.time({
        fit <- srlars::srlars(xtrain, ytrain, 
                              n_models = 1, 
                              tolerance = tolerance, 
                              x_preprocess = "ddc", 
                              y_preprocess = "wrap", 
                              cor_estimator = "wrap", 
                              cv_preprocess = "global", 
                              cv_fit = "huber", 
                              cv_loss = "huber", 
                              compute_coef = TRUE)
        preds <- predict(fit, newx = xtestdata) 
      })["elapsed"]
      
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- coef(fit)[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["RLARS",, i] <- rlars_final

    # __________________________
    # 8. Proposed: FSCRE (K=10)
    # __________________________

    fscre_final <- tryCatch({
      cpu <- system.time({
        fit <- srlars::srlars(xtrain, ytrain, 
                              n_models = n_models, 
                              tolerance = tolerance, 
                              x_preprocess = "ddc", 
                              y_preprocess = "wrap", 
                              cor_estimator = "wrap", 
                              cv_preprocess = "global", 
                              cv_fit = "huber", 
                              cv_loss = "huber", 
                              compute_coef = TRUE)
        preds <- predict(fit, newx = xtestdata) 
      })["elapsed"]
      
      mspe <- mean((preds - ytestdata)^2) / sim_data$sigma^2
      coefs <- coef(fit)[-1]
      metrics <- computeRCPR(coefs, sim_data$active_ind)
      c(MSPE = mspe, RC = metrics$rc, PR = metrics$pr, CPU = unname(cpu))
    }, error = function(e) c(NA, NA, NA, NA))
    pred_output["FSCRE",, i] <- fscre_final

    # ________________------
    # Casewise Only Methods
    # ________________------
    
    if(sim_data$contamination.scenario == "casewise"){
      
      # PENSE
      pense_final <- tryCatch({
        cpu <- system.time({
          fit <- pense::adapense_cv(x = xtrain, y = ytrain, alpha = 3/4, 
            cv_k = 5, cv_repl = 1, eps = 5e-1, explore_tol = 5e-1, 
            enpy_opts = pense::enpy_options(retain_max = 5),
            cl = cluster)
          preds <- predict(fit, xtestdata)
        })["elapsed"]
        
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
        cpu <- system.time({
          fit <- robustHD::sparseLTS(x = xtrain, y = c(ytrain), lambda = lambda_grid, mode = "lambda", cluster = cluster)
          preds <- predict(fit, xtestdata)
        })["elapsed"]
        
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