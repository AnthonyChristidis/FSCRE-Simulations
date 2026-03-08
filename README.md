# srlars

**Fast and Scalable Cellwise-Robust Ensemble (FSCRE)**

[![CRAN status](https://www.r-pkg.org/badges/version/srlars)](https://cran.r-project.org/package=srlars)

The `srlars` package implements the FSCRE algorithm, a novel, multi-stage architecture designed to perform reliable variable selection and regression in high-dimensional settings plagued by cellwise contamination.

The method establishes a robust foundation using the Detect Deviating Cells (DDC) algorithm. It then partitions the predictor space by employing a competitive ensemble architecture, where a computationally efficient, correlation-based Least Angle Regression (LARS) engine proposes candidate variables and cross-validation arbitrates their assignment.

## Installation

You can install the **stable** version from [CRAN](https://cran.r-project.org/package=srlars):

```r
install.packages("srlars")
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/srlars):

```
# install.packages("devtools")
devtools::install_github("AnthonyChristidis/srlars")
```

## Usage

This example demonstrates how to use `srlars` for robust variable selection and regression on a high-dimensional dataset with artificial cellwise contamination.

```
library(srlars)
library(mvnfast)

# 1. Simulation Parameters
n <- 50
p <- 100
p.active <- 20
snr <- 3
contamination.prop <- 0.1
set.seed(0)

# 2. Data Generation (Block Correlation)
sigma.mat <- matrix(0.2, p, p)
for(group in 0:3) sigma.mat[(group*5+1):(group*5+5), (group*5+1):(group*5+5)] <- 0.8
diag(sigma.mat) <- 1

true.beta <- c(runif(p.active, 0, 5)*(-1)^rbinom(p.active, 1, 0.7), rep(0, p - p.active))
sigma <- as.numeric(sqrt(t(true.beta) %*% sigma.mat %*% true.beta)/sqrt(snr))

x <- rmvn(n, mu = rep(0, p), sigma = sigma.mat)
y <- x %*% true.beta + rnorm(n, 0, sigma)

x_test <- rmvn(200, mu = rep(0, p), sigma = sigma.mat)
y_test <- x_test %*% true.beta + rnorm(200, 0, sigma)

# 3. Introduce Cellwise Contamination
contamination_indices <- sample(1:(n * p), round(n * p * contamination.prop))
x_train <- x
x_train[contamination_indices] <- runif(length(contamination_indices), -10, 10)

# 4. Fit the FSCRE Ensemble Model
fit <- srlars(x_train, y,
              n_models = 5,
              tolerance = 0.01,
              robust = TRUE,
              compute_coef = TRUE)

# View the disjoint sets of selected variables
print(fit$active.sets)

# 5. Prediction and Evaluation
preds <- predict(fit, newx = x_test)
mspe <- mean((y_test - preds)^2) / sigma^2
print(paste("MSPE:", round(mspe, 3)))
```

## License

This package is free and open source software, licensed under GPL (>= 2).

