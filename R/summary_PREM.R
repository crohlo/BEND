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
  class(out) <- c("summary.PREM", class(out))
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
