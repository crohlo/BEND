#' Summarize the results of a piecewise random effects model (PREM)
#'
#' @description
#' Provides a summary of a PREM model, as returned by `Bayes_PREM()`.
#'
#' @param object An object of class "PREM" (returned by `Bayes_PREM(...)`).
#' @param ... Additional arguments.
#'
#' @returns Returns a list of key parameter estimates.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # result summary
#' summary(results_prem)
#'
#' @rdname summary.PREM
#' @export
summary.PREM <- function(object, ...){

  x <- object

  # determine number of classes
  n_class <- length(unique(x$Class_Information$class_membership))

  # determine number of changepoints in each class (based on final model results)
  changepoints <- c()
  class_data <- data.frame()
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    changepoints[i] <- which.max(x$Parameter_Estimates[[class_num]]$K_prob)-1
  }
  max_cp <- max(changepoints)

  # detect covariates
  n_cov_op <- length(x$Parameter_Estimates$outcome_predictive_covariates)
  n_cov_cp <- length(x$Parameter_Estimates$class_predictive_covariates)

  # get rownames
  my_rownames <- rep(NA, 4+4*max_cp)
  my_rownames[1] = "Intercept Mean"
  my_rownames[2] = "Slope Mean"
  my_rownames[max_cp*2+2+1] = "Intercept Var"
  my_rownames[max_cp*2+2+2] = "Slope Var"
  for(k in 1:max_cp){
    my_rownames[2*k+1] = paste0("Changepoint ", k, " Mean")
    my_rownames[2*k+2] = paste0("Change in Slope ", k, " Mean")
    my_rownames[k*2+(max_cp*2+2+1)] = paste0("Changepoint ", k, " Var")
    my_rownames[k*2+(max_cp*2+2+2)] = paste0("Change in Slope ", k, " Var")
  }

  # CLASS DEPENDENT PARAMETERS -----

  # pull parameter estimates for each class
  coeff_mat_list <- vector("list", n_class)
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    cp_num <- paste0("K_", changepoints)[i]
    cur_mat <- matrix(NA, nrow=4+4*max_cp, ncol=3)
    # fixed effects
    cur_mat[1,1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[1]
    cur_mat[1,2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean_CI[1,1]
    cur_mat[1,3] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean_CI[2,1]
    cur_mat[2,1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[2]
    cur_mat[2,2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean_CI[1,2]
    cur_mat[2,3] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean_CI[2,2]
    # random effects
    cur_mat[max_cp*2+2+1,1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[1]
    cur_mat[max_cp*2+2+1,2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var_CI[1,1]
    cur_mat[max_cp*2+2+1,3] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var_CI[2,1]
    cur_mat[max_cp*2+2+2,1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[2]
    cur_mat[max_cp*2+2+2,2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var_CI[1,2]
    cur_mat[max_cp*2+2+2,3] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var_CI[2,2]
    # if there are changepoints (beyond linear)
    if(changepoints[i]>0){
      for(k in 1:changepoints[i]){
        # fixed effects
        cur_mat[2*k+1,1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_mean[k]
        cur_mat[2*k+1,2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_mean_CI[1,k]
        cur_mat[2*k+1,3] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_mean_CI[2,k]
        cur_mat[2*k+2,1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[k+2]
        cur_mat[2*k+2,2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean_CI[1,k+2]
        cur_mat[2*k+2,3] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean_CI[2,k+2]
        # random effects
        cur_mat[k*2+(max_cp*2+2+1),1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_var[k]
        cur_mat[k*2+(max_cp*2+2+1),2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_var_CI[1,k]
        cur_mat[k*2+(max_cp*2+2+1),3] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_var_CI[2,k]
        cur_mat[k*2+(max_cp*2+2+2),1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[k+2]
        cur_mat[k*2+(max_cp*2+2+2),2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var_CI[1,k+2]
        cur_mat[k*2+(max_cp*2+2+2),3] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var_CI[2,k+2]
      }
    }
    colnames(cur_mat) <- c("Estimate", "Lower CI", "Upper CI")
    rownames(cur_mat) <- paste0("Class ", i, ": ", my_rownames)
    coeff_mat_list[[i]] <- cur_mat
  }
  coeff_mat <- do.call(rbind, coeff_mat_list)

  # CLASS INDEPENDENT PARAMETERS -----
  coeff_addit <- cbind(x$Parameter_Estimates$error_var, t(x$Parameter_Estimates$error_var_CI))
  rownames(coeff_addit) <- "Error Var"

  # Covariates
  if(n_cov_op > 0){
    coeff_cov_op <- cbind(as.matrix(x$Parameter_Estimates$outcome_predictive_covariates), t(x$Parameter_Estimates$outcome_predictive_covariates_CI))
    coeff_addit <- rbind(coeff_addit, coeff_cov_op)
  }
  if(n_cov_cp > 0){
    coeff_cov_cp <- cbind(as.matrix(x$Parameter_Estimates$class_predictive_covariates), t(x$Parameter_Estimates$class_predictive_covariates_CI))
    rownames(coeff_cov_cp) <- paste0(rownames(coeff_cov_cp), " (in log-odds units)")
    coeff_log_int <- cbind(x$Parameter_Estimates$logistic_intercept, t(x$Parameter_Estimates$logistic_intercept_CI[,1]))
    rownames(coeff_log_int) <- "Logistic Intercept"
    coeff_addit <- rbind(coeff_addit, coeff_cov_cp, coeff_log_int)
  }
  colnames(coeff_addit) <- c("Estimate", "Lower CI", "Upper CI")

  ## return value
  out <- list("class_dep_params" = coeff_mat,
              "class_ind_params" = coeff_addit,
              "msrf" = x$Convergence$multivariate_psrf,
              "mean_psrf" = x$Convergence$mean_psrf,
              "DIC" = x$Model_Fit$dic)
  class(out) <- c("summary.PREM", "summary.BEND")
  return(out)
}

#' @rdname summary.PREM
#' @param x An object of class "summary.PREM".
#' @export
print.summary.PREM <- function(x, ...){

  cat("Class Dependent Parameters:\n")
  print(round(x$class_dep_params,3), na.print="")
  cat("\n")
  cat("Class Independent Parameters:\n")
  print(round(x$class_ind_params,3))
  cat("\n")
  cat("Gelman's msrf:", round(x$msrf, 3), "\n")
  cat("Mean psrf:", round(x$mean_psrf, 3), "\n")
  cat("DIC:", round(x$DIC,3))
  cat("\n")

  invisible(x)
}

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
#' @rdname summary.BPREM
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
  ss <- sub("cp", "g", gsub("cov_b_", "", rownames(ran_eff_cov_mat))) ## replace cp with g to make following code simpler
  rownames(ran_eff_cov_mat) <- paste0("Cov(Outcome ", substr(ss,1,1), ": ", param_names[substr(ss,2,2)],
                                      ", Outcome ", substr(ss,4,4), ": ", param_names[substr(ss,5,5)], ")")

  ## CORRELATION MATRIX -----
  ran_eff_corr_mat <- x$Parameter_Estimates[(2*(n_param*2)+n_covar+3+(1:n_covar)),]
  names(ran_eff_corr_mat) <- c("Estimate", "Lower CI", "Upper CI")
  # get rownames
  ss <- sub("cp", "g", gsub("corr_b_", "", rownames(ran_eff_corr_mat))) ## replace cp with g to make following code simpler
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
              "y2_var" = x$Call$y2_var,
              "fix_eff_est" = fix_eff_est,
              "ran_eff_var_mat" = ran_eff_var_mat,
              "ran_eff_cov_mat" = ran_eff_cov_mat,
              "ran_eff_corr_mat" = ran_eff_corr_mat,
              "error_cov_mat" = error_cov_mat,
              "error_corr" = error_cor_mat,
              "msrf" = x$Convergence$multivariate_psrf,
              "mean_psrf" = x$Convergence$mean_psrf,
              "DIC" = x$Model_Fit$dic)
  class(out) <- c("summary.BPREM", "summary.BEND")
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
#' @rdname summary.CREM
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
              "form" = x$Call$form,
              "fix_eff_est" = fix_eff_est,
              "ran_eff_b_var_mat" = ran_eff_b_var_mat,
              "ran_eff_b_cov_mat" = ran_eff_b_cov_mat,
              "ran_eff_g_var_mat" = ran_eff_g_var_mat,
              "ran_eff_g_cov_mat" = ran_eff_g_cov_mat,
              "error_var" = error_mat,
              "msrf" = x$Convergence$multivariate_psrf,
              "mean_psrf" = x$Convergence$mean_psrf,
              "DIC" = x$Model_Fit$dic)
  class(out) <- c("summary.CREM", "summary.BEND")
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
  cat("Functional Form:", x$form, "\n")
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
