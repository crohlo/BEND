#' Summarize the results of a piecewise random effects model (PREM)
#'
#' @description
#' Provides a summary of a PREM model, as returned by `Bayes_PREM()`.
#'
#' @param object An object of class "PREM" (returned by `Bayes_PREM(...)`).
#' @param ... Additional arguments.
#'
#' @returns Prints estimates for key parameters in the PREM. Also returns a list of these values.
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
  # Setup ----
  # determine number of classes
  n_class <- length(unique(object$Class_Information$class_membership))

  # determine number of changepoints in each class (based on final model results)
  changepoints <- c()
  class_data <- data.frame()
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    changepoints[i] <- which.max(object$Parameter_Estimates[[class_num]]$K_prob)-1
  }
  max_cp <- max(changepoints)

  # detect covariates
  n_cov_op <- length(object$Parameter_Estimates$outcome_predictive_covariates)
  n_cov_cp <- length(object$Parameter_Estimates$class_predictive_covariates)

  # CLASS DEPENDENT PARAMETERS -----
  # pull parameter estimates for each class
  coeff_mat <- matrix(NA, nrow=4+4*max_cp, ncol=n_class)
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    cp_num <- paste0("K_", changepoints)[i]
    # fixed effects
    coeff_mat[1,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[1]
    coeff_mat[2,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[2]
    # random effects
    coeff_mat[max_cp*2+2+1,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[1]
    coeff_mat[max_cp*2+2+2,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[2]
    # if there are changepoints (beyond linear)
    if(changepoints[i]>0){
      for(k in 1:changepoints[i]){
        # fixed effects
        coeff_mat[2*k+1,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_mean[k]
        coeff_mat[2*k+2,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[k+2]
        # random effects
        coeff_mat[k*2+(max_cp*2+2+1),i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_var[k]
        coeff_mat[k*2+(max_cp*2+2+2),i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_var[k+2]
      }
    }
  }
  colnames(coeff_mat) <- paste0("Class ", 1:n_class)

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
  rownames(coeff_mat) <- my_rownames

  # number of changepoints info
  coeff_n_cp <- matrix(changepoints, nrow = 1)
  rownames(coeff_n_cp) <- "Number of Changepoints"
  coeff_mat <- rbind(coeff_n_cp, coeff_mat)

  # class probability info
  class_probs <- matrix(prop.table(table(object$Class_Information$class_membership)), nrow=1)
  rownames(class_probs) <- "Empirical Class Probabilities"
  coeff_mat <- rbind(class_probs, coeff_mat)

  # CLASS INDEPENDENT PARAMETERS -----
  coeff_addit <- data.frame("Estimate" = object$Parameter_Estimates$error_var)
  rownames(coeff_addit) <- "Error Var"

  # Covariates
  if(n_cov_op > 0){
    coeff_cov_op <- data.frame("Estimate" = object$Parameter_Estimates$outcome_predictive_covariates)
    coeff_addit <- rbind(coeff_addit, coeff_cov_op)
  }
  if(n_cov_cp > 0){
    coeff_cov_cp <- data.frame("Estimate" = object$Parameter_Estimates$class_predictive_covariates)
    rownames(coeff_cov_cp) <- paste0(rownames(coeff_cov_cp), " (in log-odds units)")
    coeff_log_int <- data.frame("Estimate" = object$Parameter_Estimates$logistic_intercept)
    rownames(coeff_log_int) <- "Logistic Intercept"
    coeff_addit <- rbind(coeff_addit, coeff_cov_cp, coeff_log_int)
  }

  # PRINT output -----
  cat("Class Dependent Parameters:\n")
  print(coeff_mat, digits=3, na.print="")
  cat("\n")
  cat("Class Independent Parameters:\n")
  print(coeff_addit, digits=3)
  cat("\n")
  cat("Gelman's msrf:", round(object$Convergence$multivariate_psrf, 3), "\n")
  cat("Mean psrf:", round(object$Convergence$mean_psrf, 3), "\n")
  cat("DIC:", object$Model_Fit$dic)
  return(invisible(list("class_dep_params" = coeff_mat,
                        "class_ind_params" = coeff_addit,
                        "msrf" = object$Convergence$multivariate_psrf,
                        "mean_psrf" = object$Convergence$mean_psrf,
                        "DIC" = object$Model_Fit$dic)))
}
