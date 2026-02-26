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

  x <- object

  # determine number of parameters
  n_param <- 4

  # determine number of covariances
  n_covar <- ((n_param*2)*((n_param*2)-1))/2

  # define parameters for labeling
  param_names <- c("Intercept", "Slope", "Change in Slope", "Changepoint")
  names(param_names) <- c("0", "1", "2", "g")

  # FIXED EFFECTS -----
  fix_eff_est <- x$Parameter_Estimates[1:(n_param*2),]
  names(fix_eff_est) <- c("Estimate", "Lower CI", "Upper CI")
  rownames(fix_eff_est) <- paste0(rep(c("Outcome 1: ", "Outcome 2: "), e=4), rep(param_names,2), " Mean")

  # RANDOM EFFECTS -----

  ## COVARIANCE MATRIX -----
  ran_eff_var_mat <- x$Parameter_Estimates[((n_param*2)+(1:(n_param*2))),]
  names(ran_eff_var_mat) <- c("Estimate", "Lower CI", "Upper CI")
  rownames(ran_eff_var_mat) <- paste0(rep(c("Outcome 1: ", "Outcome 2: "), e=4), rep(param_names,2), " Var")

  ran_eff_cov_mat <- x$Parameter_Estimates[(2*(n_param*2)+(1:n_covar)),]
  names(ran_eff_cov_mat) <- c("Estimate", "Lower CI", "Upper CI")
  # get rownames
  ss <- sub("cp", "g", sub("cov_b_", "", rownames(ran_eff_cov_mat))) ## replace cp with g to make following code simpler
  rownames(ran_eff_cov_mat) <- paste0("Cov(Outcome ", substr(ss,1,1), ": ", param_names[substr(ss,2,2)],
                                      ", Outcome ", substr(ss,4,4), ": ", param_names[substr(ss,5,5)], ")")

  ## CORRELATION MATRIX -----
  ran_eff_corr_mat <- x$Parameter_Estimates[(2*(n_param*2)+n_covar+3+(1:n_covar)),]
  names(ran_eff_corr_mat) <- c("Estimate", "Lower CI", "Upper CI")
  # get rownames
  ss <- sub("cp", "g", sub("corr_b_", "", rownames(ran_eff_corr_mat))) ## replace cp with g to make following code simpler
  rownames(ran_eff_corr_mat) <- paste0("Corr(Outcome ", substr(ss,1,1), ": ", param_names[substr(ss,2,2)],
                                      ", Outcome ", substr(ss,4,4), ": ", param_names[substr(ss,5,5)], ")")

  # ERROR -----
  error_cov_mat <- x$Parameter_Estimates[(2*(n_param*2)+n_covar+(1:3)),]
  names(error_cov_mat) <- c("Estimate", "Lower CI", "Upper CI")
  rownames(error_cov_mat) <- c("Outcome 1: Error Var", "Outcome 2: Error Var", "Cov(Outcome 1: Error, Outcome 2: Error)")

  error_cor_mat <- x$Parameter_Estimates[(2*(n_param*2)+2*n_covar+3+1),]
  names(error_cor_mat) <- c("Estimate", "Lower CI", "Upper CI")
  rownames(error_cor_mat) <- c("Corr(Outcome 1: Error, Outcome 2: Error)")

  ## return value
  out <- list("data" = x$Call$data,
              "y1_var" = x$Call$y1_var,
              "y2_var" = x$Call$y1_var,
              "fix_eff_est" = fix_eff_est,
              "ran_eff_var_mat" = ran_eff_var_mat,
              "ran_eff_cov_mat" = ran_eff_cov_mat,
              "ran_eff_corr_mat" = ran_eff_corr_mat,
              "error_cov_mat" = error_cov_mat,
              "error_corr" = error_cor_mat,
              "msrf" = x$Convergence$multivariate_psrf,
              "mean_psrf" = x$Convergence$mean_psrf,
              "DIC" = x$Model_Fit$dic)
  class(out) <- c("summary.BPREM", class(out))
  return(out)
}

#' @rdname summary.BPREM
#' @param x An object of class "summary.BPREM".
#' @export
print.summary.BPREM <- function(x, ...){

  cat("Bayesian bivariate piecewise random effects model\n")
  cat("Data:", x$data, "\n")
  cat("Outcomes:", c(x$y1_var, x$y2_var), "\n")
  cat("\n")
  cat("Fixed Effect Parameters:\n")
  print(round(x$fix_eff_est,3))
  cat("\n")
  cat("Random Effect Parameters:\n")
  cat("Variances:\n")
  print(round(x$ran_eff_var_mat,3), na.print="")
  cat("\n")
  cat("Covariances:\n")
  print(round(x$ran_eff_cov_mat,3), na.print="")
  cat("\n")
  cat("Correlations:\n")
  print(round(x$ran_eff_corr_mat,3), na.print="")
  cat("\n")
  cat("Error:\n")
  cat("Variance-Covariance:\n")
  print(round(x$error_cov_mat,3), na.print="")
  cat("\n")
  cat("Correlation:\n")
  print(round(x$error_corr,3), na.print="")
  cat("\n")
  cat("Gelman's msrf:", round(x$msrf, 3), "\n")
  cat("Mean psrf:", round(x$mean_psrf, 3), "\n")
  cat("DIC:", round(x$DIC,3))
  cat("\n")

  invisible(x)
}
