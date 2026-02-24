#' Summarize the results of a bivariate piecewise random effects model (BPREM)
#'
#' @description
#' Provides a summary of a BPREM model, as returned by `Bayes_BPREM()`.
#'
#' @param object An object of class "BPREM" (returned by `Bayes_BPREM(...)`).
#' @param ... Additional arguments.
#'
#' @returns Returns a list of key parameter estimates.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_bprem)
#' # result summary
#' summary(results_bprem)
#'
#' @export
summary.BPREM <- function(object, ...){
  # Setup ----
  # determine number of parameters
  n_param <- 4

  # determine number of covariances
  n_covar <- ((n_param*2)*((n_param*2)-1))/2

  # define parameters for labeling
  param_names <- c("Intercept", "Slope", "Change in Slope", "Changepoint")

  # FIXED EFFECTS -----
  fix_eff_est <- matrix(object$Parameter_Estimates[1:(n_param*2),"Mean"], ncol=2)
  colnames(fix_eff_est) <- c("Outcome 1", "Outcome 2")
  rownames(fix_eff_est) <- paste0(param_names, " Mean")

  # RANDOM EFFECTS -----

  ## COVARIANCE MATRIX -----
  ran_eff_cov_mat <- diag(n_param*2)
  diag(ran_eff_cov_mat) <- object$Parameter_Estimates[((n_param*2)+(1:(n_param*2))),"Mean"]

  ran_eff_cov_mat[upper.tri(ran_eff_cov_mat)] <- object$Parameter_Estimates[(2*(n_param*2)+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_cov_mat[lower.tri(ran_eff_cov_mat)] <- t(ran_eff_cov_mat)[lower.tri(ran_eff_cov_mat)] # copy to lower triangle
  ran_eff_cov_mat[upper.tri(ran_eff_cov_mat)] <- NA

  colnames(ran_eff_cov_mat) <- rownames(ran_eff_cov_mat) <- paste0(rep(c("Outcome 1: ", "Outcome 2: "), e=4), rep(param_names,2))

  ## CORRELATION MATRIX -----
  ran_eff_corr_mat <- diag(n_param*2)

  ran_eff_corr_mat[upper.tri(ran_eff_corr_mat)] <- object$Parameter_Estimates[(2*(n_param*2)+n_covar+3+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_corr_mat[lower.tri(ran_eff_corr_mat)] <- t(ran_eff_corr_mat)[lower.tri(ran_eff_corr_mat)] # copy to lower triangle
  ran_eff_corr_mat[upper.tri(ran_eff_corr_mat)] <- NA

  colnames(ran_eff_corr_mat) <- rownames(ran_eff_corr_mat) <- paste0(rep(c("Outcome 1: ", "Outcome 2: "), e=4), rep(param_names,2))

  # ERROR -----
  error_cov_mat <- diag(2)
  diag(error_cov_mat) <- object$Parameter_Estimates[(2*(n_param*2)+n_covar+(1:2)),"Mean"]

  error_cov_mat[upper.tri(error_cov_mat)] <- error_cov_mat[lower.tri(error_cov_mat)] <- object$Parameter_Estimates[(2*(n_param*2)+n_covar+(3)),"Mean"] # covariances on upper triangle
  error_cov_mat[upper.tri(error_cov_mat)] <- NA

  colnames(error_cov_mat) <- rownames(error_cov_mat) <- paste0(c("Outcome 1: ", "Outcome 2: "), "Error")

  ## return value
  out <- list("fix_eff_est" = fix_eff_est,
              "ran_eff_cov_mat" = ran_eff_cov_mat,
              "ran_eff_corr_mat" = ran_eff_corr_mat,
              "error_cov_mat" = error_cov_mat,
              "error_corr" = object$Parameter_Estimates[(2*(n_param*2)+2*n_covar+3+1),"Mean"],
              "msrf" = object$Convergence$multivariate_psrf,
              "mean_psrf" = object$Convergence$mean_psrf,
              "DIC" = object$Model_Fit$dic)
  class(out) <- c("summary.BPREM", class(object))
  return(out)
}

#' @rdname summary.BPREM
#' @export
print.summary.BPREM <- function(object, ...){

  cat("Fixed Effect Parameters:\n")
  print(round(object$fix_eff_est,3))
  cat("\n")
  cat("Random Effect Parameters:\n")
  cat("Covariance Matrix:\n")
  print(round(object$ran_eff_cov_mat,3), na.print="")
  cat("\n")
  cat("Correlation Matrix:\n")
  print(round(object$ran_eff_corr_mat,3), na.print="")
  cat("\n")
  cat("Error:\n")
  cat("Covariance Matrix:\n")
  print(round(object$error_cov_mat,3), na.print="")
  cat("\n")
  cat("Error Corr:", round(object$error_corr,3), "\n")
  cat("Gelman's msrf:", round(object$msrf, 3), "\n")
  cat("Mean psrf:", round(object$mean_psrf, 3), "\n")
  cat("DIC:", round(object$DIC,3))
  cat("\n")

  invisible(object)
}
