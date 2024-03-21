#' Summarize the results of a crossed random effects model (CREM)
#'
#' @description
#' Provides a summary of a CREM model, as returned by `Bayes_CREM()`.
#'
#' @param object An object of class "CREM" (returned by `Bayes_CREM(...)`).
#' @param ... Additional arguments.
#'
#' @returns Prints estimates for key parameters in the CREM. Also returns a list of these values.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_pcrem)
#' # result summary
#' summary(results_pcrem)
#'
#' @export
summary.CREM <- function(object, ...){
  # Setup ----
  # determine form
  form <- object$Functional_Form

  # determine number of parameters
  if(form=="linear")                           n_param <- 2
  if(form=="quadratic" | form=="exponential")  n_param <- 3
  if(form=="piecewise")                        n_param <- 4

  # determine number of covariances
  n_covar <- (n_param*(n_param-1))/2

  # define parameters for labeling
  if(form=="linear")       param_names <- c("Intercept", "Slope")
  if(form=="quadratic")    param_names <- c("Intercept", "Linear Slope", "Quadratic Slope")
  if(form=="exponential")  param_names <- c("Intercept", "Total Change", "Growth Rate")
  if(form=="piecewise")    param_names <- c("Intercept", "Slope", "Change in Slope", "Changepoint")

  # FIXED EFFECTS -----
  fix_eff_est <- matrix(object$Parameter_Estimates[1:n_param,"Mean"], ncol=1)
  colnames(fix_eff_est) <- "Estimate"
  rownames(fix_eff_est) <- paste0(param_names, " Mean")

  # RANDOM EFFECTS -----

  ## INDIVIDUALS -----
  ran_eff_b_mat <- diag(n_param)
  diag(ran_eff_b_mat) <- object$Parameter_Estimates[(n_param+(1:n_param)),"Mean"]

  ran_eff_b_mat[upper.tri(ran_eff_b_mat)] <- object$Parameter_Estimates[(2*n_param+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_b_mat[lower.tri(ran_eff_b_mat)] <- t(ran_eff_b_mat)[lower.tri(ran_eff_b_mat)] # copy to lower triangle
  ran_eff_b_mat[upper.tri(ran_eff_b_mat)] <- NA

  colnames(ran_eff_b_mat) <- rownames(ran_eff_b_mat) <- param_names

  ## GROUPS -----
  ran_eff_g_mat <- diag(n_param)
  diag(ran_eff_g_mat) <- object$Parameter_Estimates[(2*n_param+n_covar+(1:n_param)),"Mean"]

  ran_eff_g_mat[upper.tri(ran_eff_g_mat)] <- object$Parameter_Estimates[(3*n_param+n_covar+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_g_mat[lower.tri(ran_eff_g_mat)] <- t(ran_eff_g_mat)[lower.tri(ran_eff_g_mat)] # copy to lower triangle
  ran_eff_g_mat[upper.tri(ran_eff_g_mat)] <- NA

  colnames(ran_eff_g_mat) <- rownames(ran_eff_g_mat) <- param_names

  # PRINT output -----
  cat("Fixed Effect Parameters:\n")
  print(fix_eff_est, digits=3)
  cat("\n")
  cat("Random Effect Parameters:\n")
  cat("Individual Random Effects Covariance Matrix:\n")
  print(ran_eff_b_mat, digits=3, na.print="")
  cat("\n")
  cat("Group Random Effects Covariance Matrix:\n")
  print(ran_eff_g_mat, digits=3, na.print="")
  cat("\n")
  cat("Error Var:", object$Parameter_Estimates[(3*n_param+2*n_covar+1),"Mean"], "\n")
  cat("Gelman's msrf:", round(object$Convergence$multivariate_psrf, 3), "\n")
  cat("Mean psrf:", round(object$Convergence$mean_psrf, 3), "\n")
  cat("DIC:", object$Model_Fit$dic)
  return(invisible(list("fix_eff_est" = fix_eff_est,
                        "ran_eff_b_mat" = ran_eff_b_mat,
                        "ran_eff_g_mat" = ran_eff_g_mat,
                        "error_var" = object$Parameter_Estimates[(3*n_param+2*n_covar+1),"Mean"],
                        "msrf" = object$Convergence$multivariate_psrf,
                        "mean_psrf" = object$Convergence$mean_psrf,
                        "DIC" = object$Model_Fit$dic)))
}
