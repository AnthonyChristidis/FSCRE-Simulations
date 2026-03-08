#' @description computeRCPR checks for non-zero coefficients and returns the recall and precision of an estimator.
#' @param coef_vector Vector of coefficients.
#' @param active_preds Number of active predictors.

computeRCPR <- function(coefs, active_ind) {
    active_set_hat <- which(coefs != 0)
    if (length(active_set_hat) == 0) return(list(pr = 0, rc = 0))
    pr <- length(intersect(active_set_hat, active_ind)) / length(active_set_hat)
    rc <- length(intersect(active_set_hat, active_ind)) / length(active_ind)
    return(list(pr = pr, rc = rc))
  }
