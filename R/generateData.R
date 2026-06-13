#' Generate Contaminated and Uncontaminated Data for Simulation
#' 
#' @description \code{generateData} produces contaminated training data sets and a large uncontaminated test set
#' based on a block-diagonal correlation structure and various contamination scenarios.
#' 
#' @param N Number of training sets (replications).
#' @param n Sample size for each training set.
#' @param m Sample size for the test set.
#' @param p Total number of predictors.
#' @param rho Correlation within a block of active parameters.
#' @param rho.inactive Correlation between blocks of active parameters.
#' @param p.active Number of active parameters (true signals).
#' @param group.size Size of one block of active parameters.
#' @param snr Signal-to-noise ratio.
#' @param contamination.prop Contamination proportion (scalar or vector of length 2 for mixtures).
#' @param contamination.scenario Type of contamination: "casewise", "cellwise_marginal", "cellwise_correlation", "mixture_marginal", or "mixture_correlation".
#'
#' @return A list containing training data, test data, and true model parameters.
generateData <- function(N, 
                         n, 
                         m, 
                         p, 
                         rho, 
                         rho.inactive, 
                         p.active, 
                         group.size, 
                         snr, 
                         contamination.prop,
                         contamination.scenario){ 
  
  # Lists to store training data
  xlist <- list()
  ylist <- list()
  
  # Casewise contamination parameters
  k_lev <- 2
  k_slo <- 100
  # Cellwise contamination parameters
  gamma <- 3
  cell.mean <- 10
  
  # Setting up correlation between and within blocks of active parameters
  sigma.mat <- matrix(0, nrow = p, ncol = p)
  sigma.mat[1:p.active, 1:p.active] <- rho.inactive
  for(group in 0:(p.active/group.size - 1)) {
    sigma.mat[(group*group.size+1):(group*group.size+group.size), 
              (group*group.size+1):(group*group.size+group.size)] <- rho
  }
  diag(sigma.mat) <- 1
  
  # Generating true coefficient vector (sparse)
  trueBeta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
  
  # Contamination of trueBeta for casewise scenarios
  beta_cont <- trueBeta
  beta_cont[trueBeta!=0] <- beta_cont[trueBeta!=0]*(1 + k_slo)
  beta_cont[trueBeta==0] <- k_slo*max(abs(trueBeta))
  
  # Noise parameter scaled to target SNR
  sigma <- as.numeric(sqrt(t(trueBeta) %*% sigma.mat %*% trueBeta)/sqrt(snr))
  
  # Simulating clean training data baseline
  for(i in 1:N){
    x_train <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigma.mat)
    y <- x_train %*% trueBeta + rnorm(n, 0, sigma)
    
    xlist[[i]] <- x_train
    ylist[[i]] <- y
  }
  
  # Apply Contamination Scenarios
  if(contamination.scenario == "casewise") {
    
    contamination_indices <- 1:floor(n*contamination.prop)
    for(i in 1:N) {
      for(cont_id in contamination_indices){
        # Generate leverage point direction
        a <- runif(p, min = -1, max = 1)
        a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
        xlist[[i]][cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
          k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
        # Generate vertical outlier
        ylist[[i]][cont_id] <- t(xlist[[i]][cont_id,]) %*% beta_cont
      }
    }
      
  } else if(contamination.scenario == "cellwise_marginal") {
      
    for(i in 1:N) {
      x_train <- xlist[[i]]
      # Contaminate ONLY the predictors (X), leaving y clean
      contamination_indices <- sample.int(n * p, round(n * p * contamination.prop))
      x_train[contamination_indices] <- rnorm(length(contamination_indices), cell.mean, 1)
      
      xlist[[i]] <- x_train
      # ylist[[i]] remains unchanged
    }
      
  } else if(contamination.scenario == "cellwise_correlation") {
      
    for(i in 1:N) {
      contamination_indices <- sample.int(n * p, round(n * p * contamination.prop))
      x_train <- xlist[[i]]
      x_train[contamination_indices] <- NA
      for(row_id in 1:n){
        cells_id <- which(is.na(x_train[row_id,]))
        if(length(cells_id) > 0) {
          if (length(cells_id) > 1) {
            # Multivariate correlation outlier
            mu_cells <- rep(0, length(cells_id))
            sigma_cells <- sigma.mat[cells_id, cells_id, drop=FALSE]
            # Use symmetric=TRUE for safety
            eigen_vec <- eigen(sigma_cells, symmetric = TRUE)$vectors[, length(cells_id)]
            denom <- as.numeric(mahalanobis(t(eigen_vec), center = mu_cells, cov = sigma_cells))
            
            x_train[row_id, cells_id] <- gamma * sqrt(length(cells_id)) * t(eigen_vec) / sqrt(denom)
          } else {
            # Fallback for a single cell outlier
            x_train[row_id, cells_id] <- gamma * 3 
          }
        }
      }
      xlist[[i]] <- x_train
    }
      
  } else if(contamination.scenario == "mixture_marginal"){
      
    n.casewise <- floor(n*contamination.prop[1])
    
    for(i in 1:N) {
      # Casewise Contamination (First chunk of rows)
      if (n.casewise > 0) {
        contamination_indices <- 1:n.casewise
        for(cont_id in contamination_indices){
          a <- runif(p, min = -1, max = 1)
          a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
          xlist[[i]][cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
            k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
          ylist[[i]][cont_id] <- t(xlist[[i]][cont_id,]) %*% beta_cont
        }
      }
    
      # Cellwise Marginal Contamination (Remaining rows) - X ONLY
      if (n - n.casewise > 0) {
        subset_idx <- (n.casewise + 1):n
        sub_x <- xlist[[i]][subset_idx, , drop = FALSE]
        n_sub <- nrow(sub_x)
        
        contamination_indices <- sample.int(n_sub * p, round(n_sub * p * contamination.prop[2]))
        sub_x[contamination_indices] <- rnorm(length(contamination_indices), cell.mean, 1)
        
        xlist[[i]][subset_idx, ] <- sub_x
        # ylist[[i]] remains unchanged for these rows
      }
    }
    
  } else if(contamination.scenario == "mixture_correlation"){
      
    n.casewise <- floor(n*contamination.prop[1])
    
    for(i in 1:N) {
      # Casewise Contamination
      if (n.casewise > 0) {
        contamination_indices <- 1:n.casewise
        for(cont_id in contamination_indices){
          a <- runif(p, min = -1, max = 1)
          a <- a - as.numeric((1/p)*t(a) %*% rep(1, p))
          xlist[[i]][cont_id,] <- mvnfast::rmvn(1, rep(0, p), 0.1^2*diag(p)) + 
            k_lev * a / as.numeric(sqrt(t(a) %*% solve(sigma.mat) %*% a))
          ylist[[i]][cont_id] <- t(xlist[[i]][cont_id,]) %*% beta_cont
        }
      }
      
      # Cellwise Correlation Contamination (Remaining rows)
      if (n - n.casewise > 0) {
        subset_idx <- (n.casewise + 1):n
        sub_matrix <- xlist[[i]][subset_idx, , drop = FALSE]
        n_sub <- nrow(sub_matrix)
        
        contamination_indices <- sample.int(n_sub * p, round(n_sub * p * contamination.prop[2]))
        sub_matrix[contamination_indices] <- NA
        
        for(row_id in 1:n_sub){
          cells_id <- which(is.na(sub_matrix[row_id,]))
          if(length(cells_id) > 0) {
            if (length(cells_id) > 1) {
              mu_cells <- rep(0, length(cells_id))
              sigma_cells <- sigma.mat[cells_id, cells_id, drop=FALSE]
              eigen_vec <- eigen(sigma_cells, symmetric = TRUE)$vectors[, length(cells_id)]
              denom <- as.numeric(mahalanobis(t(eigen_vec), center = mu_cells, cov = sigma_cells))
              
              sub_matrix[row_id, cells_id] <- gamma * sqrt(length(cells_id)) * t(eigen_vec) / sqrt(denom)
            } else {
              sub_matrix[row_id, cells_id] <- gamma * 3
            }
          }
        }
        xlist[[i]][subset_idx, ] <- sub_matrix
      }
    }
  }
  
  # Simulating uncontaminated test data (massive size m)
  x_test <- mvnfast::rmvn(m, mu = rep(0, p), sigma = sigma.mat)
  y_test <- x_test %*% trueBeta + rnorm(m, 0, sigma)
  
  # Ensure y vectors are numeric vectors, not matrices
  for (i in 1:N) {
      ylist[[i]] <- as.numeric(ylist[[i]])
  }
  y_test <- as.numeric(y_test)
  
  return(
   list(training_data = list(xtrain = xlist, ytrain = ylist), 
        testing_data = list(xtest = x_test, ytest = y_test), 
        pactive = p.active, n = n, sigma = sigma, 
        active_ind = which(trueBeta != 0), p = p, 
        trueBeta = as.numeric(trueBeta),
        contamination.scenario = contamination.scenario)
  )
}