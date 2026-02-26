#' Print the results of a crossed random effects model (CREM)
#'
#' @description
#' Provides a summary of a CREM model, as returned by `Bayes_CREM()`.
#'
#' @param x An object of class "CREM" (returned by `Bayes_CREM(...)`).
#' @param ... Additional arguments.
#'
#' @returns Returns a list of key parameter estimates.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_pcrem)
#' # print results
#' print(results_pcrem)
#'
#' @export
print.CREM <- function(x, ...){

  # determine form
  form <- x$Functional_Form

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
  fix_eff_est <- matrix(x$Parameter_Estimates[1:n_param,"Mean"], ncol=1)
  colnames(fix_eff_est) <- "Estimate"
  rownames(fix_eff_est) <- paste0(param_names, " Mean")

  # RANDOM EFFECTS -----

  ## INDIVIDUALS -----
  ran_eff_b_mat <- diag(n_param)
  diag(ran_eff_b_mat) <- x$Parameter_Estimates[(n_param+(1:n_param)),"Mean"]

  ran_eff_b_mat[upper.tri(ran_eff_b_mat)] <- x$Parameter_Estimates[(2*n_param+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_b_mat[lower.tri(ran_eff_b_mat)] <- t(ran_eff_b_mat)[lower.tri(ran_eff_b_mat)] # copy to lower triangle
  ran_eff_b_mat[upper.tri(ran_eff_b_mat)] <- NA

  colnames(ran_eff_b_mat) <- rownames(ran_eff_b_mat) <- param_names

  ## GROUPS -----
  ran_eff_g_mat <- diag(n_param)
  diag(ran_eff_g_mat) <- x$Parameter_Estimates[(2*n_param+n_covar+(1:n_param)),"Mean"]

  ran_eff_g_mat[upper.tri(ran_eff_g_mat)] <- x$Parameter_Estimates[(3*n_param+n_covar+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_g_mat[lower.tri(ran_eff_g_mat)] <- t(ran_eff_g_mat)[lower.tri(ran_eff_g_mat)] # copy to lower triangle
  ran_eff_g_mat[upper.tri(ran_eff_g_mat)] <- NA

  colnames(ran_eff_g_mat) <- rownames(ran_eff_g_mat) <- param_names

  ## return value
  out <- list("fix_eff_est" = fix_eff_est,
              "ran_eff_b_mat" = ran_eff_b_mat,
              "ran_eff_g_mat" = ran_eff_g_mat,
              "error_var" = x$Parameter_Estimates[(3*n_param+2*n_covar+1),"Mean"],
              "msrf" = x$Convergence$multivariate_psrf,
              "mean_psrf" = x$Convergence$mean_psrf,
              "DIC" = x$Model_Fit$dic)

  cat("Bayesian crossed random effects model\n")
  cat("Data:", x$Call$data, "\n")
  cat("Outcome:", x$Call$y_var, "\n")
  cat("Individuals:", x$Call$ind_id_var, "\n")
  cat("Group:", x$Call$cross_id_var, "\n")
  cat("\n")
  cat("Fixed Effect Parameters:\n")
  print(round(out$fix_eff_est,3))
  cat("\n")
  cat("Random Effect Parameters:\n")
  cat("Individual Random Effects Covariance Matrix:\n")
  print(round(out$ran_eff_b_mat,3), na.print="")
  cat("\n")
  cat("Group Random Effects Covariance Matrix:\n")
  print(round(out$ran_eff_g_mat,3), na.print="")
  cat("\n")
  cat("Error Var:", round(out$error_var,3), "\n")
  cat("Gelman's msrf:", round(out$msrf, 3), "\n")
  cat("Mean psrf:", round(out$mean_psrf, 3), "\n")
  cat("DIC:", round(out$DIC,3))
  cat("\n")

  invisible(x)
}

