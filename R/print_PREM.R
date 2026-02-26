#' Print the results of a piecewise random effects model (PREM)
#'
#' @description
#' Provides a summary of a PREM model, as returned by `Bayes_PREM()`.
#'
#' @param x An object of class "PREM" (returned by `Bayes_PREM(...)`).
#' @param ... Additional arguments.
#'
#' @returns Returns a list of key parameter estimates.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # print results
#' print(results_prem)
#'
#' @export
print.PREM <- function(x, ...){

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

  # CLASS DEPENDENT PARAMETERS -----
  # pull fixed effects for each class
  coeff_mat <- matrix(NA, nrow=2+2*max_cp, ncol=n_class)
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    cp_num <- paste0("K_", changepoints)[i]
    # fixed effects
    coeff_mat[1,i] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[1]
    coeff_mat[2,i] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[2]
    # if there are changepoints (beyond linear)
    if(changepoints[i]>0){
      for(k in 1:changepoints[i]){
        # fixed effects
        coeff_mat[2*k+1,i] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_mean[k]
        coeff_mat[2*k+2,i] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[k+2]
      }
    }
  }
  colnames(coeff_mat) <- paste0("Class ", 1:n_class)

  my_rownames <- rep(NA, 2+2*max_cp)
  my_rownames[1] = "Intercept"
  my_rownames[2] = "Slope"
  for(k in 1:max_cp){
    my_rownames[2*k+1] = paste0("Changepoint ", k)
    my_rownames[2*k+2] = paste0("Change in Slope ", k)
  }
  rownames(coeff_mat) <- paste0(my_rownames, " Mean")

  # pull variances for each class
  varcor_list <- vector("list", n_class)
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    cp_num <- paste0("K_", changepoints)[i]
    varcor_vec <- c()
    # random effects
    varcor_vec[1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[1]
    varcor_vec[2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[2]
    # if there are changepoints (beyond linear)
    if(changepoints[i]>0){
      for(k in 1:changepoints[i]){
        # random effects
        varcor_vec[2*k+1] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_var[k]
        varcor_vec[2*k+2] <- x$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[k+2]
      }
    }
    n_var <- length(varcor_vec)

    # param_names
    param_names <- rep(NA, n_var)
    param_names[1] = "Intercept"
    param_names[2] = "Slope"
    if(changepoints[i]>0){
      for(k in 1:changepoints[i]){
        param_names[2*k+1] = paste0("Changepoint ", k)
        param_names[2*k+2] = paste0("Change in Slope ", k)
      }
    }
    # format matrix
    varcor_mat <- diag(n_var)
    diag(varcor_mat) <- varcor_vec
    varcor_mat[row(varcor_mat) != col(varcor_mat)] <- NA
    colnames(varcor_mat) <- rownames(varcor_mat) <- param_names

    # save to list
    varcor_list[[i]] <- varcor_mat
  }
  names(varcor_list) <- paste0("class_", 1:n_class)

  # number of changepoints info
  coeff_n_cp <- matrix(changepoints, nrow = 1)
  rownames(coeff_n_cp) <- "Number of Changepoints"
  coeff_mat <- rbind(coeff_n_cp, coeff_mat)

  # class probability info
  class_probs <- matrix(prop.table(table(x$Class_Information$class_membership)), nrow=1)
  rownames(class_probs) <- "Empirical Class Probabilities"
  coeff_mat <- rbind(class_probs, coeff_mat)

  # CLASS INDEPENDENT PARAMETERS -----
  coeff_addit <- data.frame("Estimate" = x$Parameter_Estimates$error_var)
  rownames(coeff_addit) <- "Error Var"

  # Covariates
  if(n_cov_op > 0){
    coeff_cov_op <- data.frame("Estimate" = x$Parameter_Estimates$outcome_predictive_covariates)
    coeff_addit <- rbind(coeff_addit, coeff_cov_op)
  }
  if(n_cov_cp > 0){
    coeff_cov_cp <- data.frame("Estimate" = x$Parameter_Estimates$class_predictive_covariates)
    rownames(coeff_cov_cp) <- paste0(rownames(coeff_cov_cp), " (in log-odds units)")
    coeff_log_int <- data.frame("Estimate" = x$Parameter_Estimates$logistic_intercept)
    rownames(coeff_log_int) <- "Logistic Intercept"
    coeff_addit <- rbind(coeff_addit, coeff_cov_cp, coeff_log_int)
  }

  ## return value
  out <- list("class_dep_params" = coeff_mat,
              "class_dep_varcov" = varcor_list,
              "class_ind_params" = coeff_addit,
              "msrf" = x$Convergence$multivariate_psrf,
              "mean_psrf" = x$Convergence$mean_psrf,
              "DIC" = x$Model_Fit$dic)

  cat("Bayesian piecewise random effects model\n")
  cat("Data:", x$Call$data, "\n")
  cat("Outcome:", x$Call$y_var, "\n")
  cat("Number of Classes:", n_class, "\n")
  cat("Maximum Number of Changepoints:", n_class, "\n")
  if(n_cov_cp==1) cat("Class-Predictive Covariates:", x$Call$class_predictive_vars,"\n")
  if(n_cov_cp>1) cat("Class-Predictive Covariates:", as.character(x$Call$class_predictive_vars)[-1],"\n")
  if(n_cov_op==1) cat("Outcome-Predictive Covariates:", x$Call$outcome_predictive_vars,"\n")
  if(n_cov_op>1) cat("Outcome-Predictive Covariates:", as.character(x$Call$outcome_predictive_vars)[-1],"\n")
  cat("\n")
  cat("Class Dependent Parameters:\n")
  cat("Fixed Effects Parameters:\n")
  print(round(out$class_dep_params,3), na.print="")
  cat("\n")
  cat("Random Effects Parameters:\n")
  for(i in 1:length(out$class_dep_varcov)){
    cat(paste0("Class ",i,":\n"))
    print(round(out$class_dep_varcov[[i]],3), na.print="")
    cat("\n")
  }
  cat("Class Independent Parameters:\n")
  print(round(out$class_ind_params,3))
  cat("\n")
  cat("Gelman's msrf:", round(out$msrf, 3), "\n")
  cat("Mean psrf:", round(out$mean_psrf, 3), "\n")
  cat("DIC:", round(out$DIC,3))
  cat("\n")

  invisible(x)
}
