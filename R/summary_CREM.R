#' Summarize the results of a crossed random effects model (CREM)
#'
#' @description
#' Provides a summary of a CREM model, as returned by `Bayes_CREM()`.
#'
#' @param object An object of class "CREM" (returned by `Bayes_CREM(...)`).
#' @param ... Additional arguments.
#'
#' @returns Returns a list of key parameter estimates.
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

  x <- object

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
  if(form=="linear")       names(param_names) <- c("0", "1")
  if(form=="quadratic")    param_names <- c("Intercept", "Linear Slope", "Quadratic Slope")
  if(form=="quadratic")    names(param_names) <- c("0", "1", "2")
  if(form=="exponential")  param_names <- c("Intercept", "Total Change", "Growth Rate")
  if(form=="exponential")  names(param_names) <- c("0", "1", "2")
  if(form=="piecewise")    param_names <- c("Intercept", "Slope", "Change in Slope", "Changepoint")
  if(form=="piecewise")    names(param_names) <- c("0", "1", "2", "g")

  # FIXED EFFECTS -----
  fix_eff_est <- x$Parameter_Estimates[1:n_param,]
  colnames(fix_eff_est) <- c("Estimate", "Lower CI", "Upper CI")
  rownames(fix_eff_est) <- paste0(param_names, " Mean")

  # RANDOM EFFECTS -----

  ## INDIVIDUALS -----
  ran_eff_b_var_mat <- x$Parameter_Estimates[(n_param+(1:n_param)),]
  colnames(ran_eff_b_var_mat) <- c("Estimate", "Lower CI", "Upper CI")
  rownames(ran_eff_b_var_mat) <- paste0(param_names, " Var")

  ran_eff_b_cov_mat <- x$Parameter_Estimates[(2*n_param+(1:n_covar)),]
  colnames(ran_eff_b_cov_mat) <- c("Estimate", "Lower CI", "Upper CI")
  # get rownames
  ss <- sub("cp", "g", sub("cov_b_", "", rownames(ran_eff_b_cov_mat))) ## replace cp with g to make following code simpler
  rownames(ran_eff_b_cov_mat) <- paste0("Cov(", param_names[substr(ss,1,1)], ", ", param_names[substr(ss,2,2)], ")")

  ## GROUPS -----
  ran_eff_g_var_mat <- x$Parameter_Estimates[(2*n_param+n_covar+(1:n_param)),]
  colnames(ran_eff_g_var_mat) <- c("Estimate", "Lower CI", "Upper CI")
  rownames(ran_eff_g_var_mat) <- paste0(param_names, " Var")

  ran_eff_g_cov_mat <- x$Parameter_Estimates[(3*n_param+n_covar+(1:n_covar)),]
  colnames(ran_eff_g_cov_mat) <- c("Estimate", "Lower CI", "Upper CI")
  # get rownames
  ss <- sub("cp", "g", sub("cov_g_", "", rownames(ran_eff_g_cov_mat))) ## replace cp with g to make following code simpler
  rownames(ran_eff_g_cov_mat) <- paste0("Cov(", param_names[substr(ss,1,1)], ", ", param_names[substr(ss,2,2)], ")")

  ## ERROR -----
  error_mat <- x$Parameter_Estimates[(3*n_param+2*n_covar+1),]
  colnames(error_mat) <- c("Estimate", "Lower CI", "Upper CI")
  rownames(error_mat) <- c("Error Var")

  ## return value
  out <- list("data" = x$Call$data,
              "y_var" = x$Call$y_var,
              "ind_id_var" = x$Call$ind_id_var,
              "cross_id_var" = x$Call$cross_id_var,
              "fix_eff_est" = fix_eff_est,
              "ran_eff_b_var_mat" = ran_eff_b_var_mat,
              "ran_eff_b_cov_mat" = ran_eff_b_cov_mat,
              "ran_eff_g_var_mat" = ran_eff_g_var_mat,
              "ran_eff_g_cov_mat" = ran_eff_g_cov_mat,
              "error_var" = error_mat,
              "msrf" = x$Convergence$multivariate_psrf,
              "mean_psrf" = x$Convergence$mean_psrf,
              "DIC" = x$Model_Fit$dic)
  class(out) <- c("summary.CREM", class(out))
  return(out)
}

#' @rdname summary.CREM
#' @param x An object of class "summary.CREM".
#' @export
print.summary.CREM <- function(x, ...){

  cat("Bayesian crossed random effects model\n")
  cat("Data:", x$data, "\n")
  cat("Outcome:", x$y_var, "\n")
  cat("Individuals:", x$ind_id_var, "\n")
  cat("Group:", x$cross_id_var, "\n")
  cat("\n")
  cat("Fixed Effect Parameters:\n")
  print(round(x$fix_eff_est,3))
  cat("\n")
  cat("Random Effect Parameters:\n")
  cat("Individual Random Effects Variance-Covariance:\n")
  print(round(rbind(x$ran_eff_b_var_mat, x$ran_eff_b_cov_mat),3), na.print="")
  cat("\n")
  cat("Group Random Effects Variance-Covariance:\n")
  print(round(rbind(x$ran_eff_g_var_mat, x$ran_eff_g_cov_mat),3), na.print="")
  cat("\n")
  cat("Error Variance:\n")
  print(round(x$error_var,3), na.print="")
  cat("\n")
  cat("Gelman's msrf:", round(x$msrf, 3), "\n")
  cat("Mean psrf:", round(x$mean_psrf, 3), "\n")
  cat("DIC:", round(x$DIC,3))
  cat("\n")

  invisible(x)
}
