#' Extract random effects variance-covariance matrix
#'
#' @description
#' Extracts the random effects variance-covariance matrix from a fitted model of class "BPREM", "CREM", or "PREM".
#'
#' @param x An object of class "BPREM", "CREM", or "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a list of the random effects variance-covariance matrices.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get random effects variance-covariance matrices
#' getVarCov(results_prem)
#'
#' @export
getVarCov <- function(x, ...) UseMethod(("getVarCov"))

#' @rdname getVarCov
#' @export
getVarCov.BPREM <- function(x, ...){

  # determine number of parameters
  n_param <- 4
  # determine number of covariances
  n_covar <- ((n_param*2)*((n_param*2)-1))/2
  # define parameters for labeling
  param_names <- c("Intercept", "Slope", "Change in Slope", "Changepoint")

  ## COVARIANCE MATRIX
  ran_eff_cov_mat <- diag(n_param*2)
  diag(ran_eff_cov_mat) <- x$Parameter_Estimates[((n_param*2)+(1:(n_param*2))),"Mean"]

  ran_eff_cov_mat[upper.tri(ran_eff_cov_mat)] <- x$Parameter_Estimates[(2*(n_param*2)+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_cov_mat[lower.tri(ran_eff_cov_mat)] <- t(ran_eff_cov_mat)[lower.tri(ran_eff_cov_mat)] # copy to lower triangle
  ran_eff_cov_mat[upper.tri(ran_eff_cov_mat)] <- NA

  colnames(ran_eff_cov_mat) <- rownames(ran_eff_cov_mat) <- paste0(rep(paste0(c(x$Call$y1_var, x$Call$y2_var), ": "), e=4), rep(param_names,2))

  ## CORRELATION MATRIX
  ran_eff_corr_mat <- diag(n_param*2)

  ran_eff_corr_mat[upper.tri(ran_eff_corr_mat)] <- x$Parameter_Estimates[(2*(n_param*2)+n_covar+3+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_corr_mat[lower.tri(ran_eff_corr_mat)] <- t(ran_eff_corr_mat)[lower.tri(ran_eff_corr_mat)] # copy to lower triangle
  ran_eff_corr_mat[upper.tri(ran_eff_corr_mat)] <- NA

  colnames(ran_eff_corr_mat) <- rownames(ran_eff_corr_mat) <- paste0(rep(paste0(c(x$Call$y1_var, x$Call$y2_var), ": "), e=4), rep(param_names,2))

  out <- list(cov_mat = ran_eff_cov_mat,
              cor_mat = ran_eff_corr_mat)
  class(out) <- c("getVarCov.BPREM", class(out))
  return(out)
}

#' @rdname getVarCov
#' @export
getVarCov.CREM <- function(x, ...){

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

  ## INDIVIDUALS
  ran_eff_b_mat <- diag(n_param)
  diag(ran_eff_b_mat) <- x$Parameter_Estimates[(n_param+(1:n_param)),"Mean"]

  ran_eff_b_mat[upper.tri(ran_eff_b_mat)] <- x$Parameter_Estimates[(2*n_param+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_b_mat[lower.tri(ran_eff_b_mat)] <- t(ran_eff_b_mat)[lower.tri(ran_eff_b_mat)] # copy to lower triangle
  ran_eff_b_mat[upper.tri(ran_eff_b_mat)] <- NA

  colnames(ran_eff_b_mat) <- rownames(ran_eff_b_mat) <- param_names

  ## GROUPS
  ran_eff_g_mat <- diag(n_param)
  diag(ran_eff_g_mat) <- x$Parameter_Estimates[(2*n_param+n_covar+(1:n_param)),"Mean"]

  ran_eff_g_mat[upper.tri(ran_eff_g_mat)] <- x$Parameter_Estimates[(3*n_param+n_covar+(1:n_covar)),"Mean"] # covariances on upper triangle
  ran_eff_g_mat[lower.tri(ran_eff_g_mat)] <- t(ran_eff_g_mat)[lower.tri(ran_eff_g_mat)] # copy to lower triangle
  ran_eff_g_mat[upper.tri(ran_eff_g_mat)] <- NA

  colnames(ran_eff_g_mat) <- rownames(ran_eff_g_mat) <- param_names

  out <- list(ind_mat = ran_eff_b_mat,
              grp_mat = ran_eff_g_mat)
  class(out) <- c("getVarCov.CREM", class(out))
  return(out)
}

#' @rdname getVarCov
#' @export
getVarCov.PREM <- function(x, ...){

  # determine number of classes
  n_class <- length(unique(x$Class_Information$class_membership))

  # determine number of changepoints in each class (based on final model results)
  changepoints <- c()
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    changepoints[i] <- which.max(x$Parameter_Estimates[[class_num]]$K_prob)-1
  }

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

  out <- varcor_list
  class(out) <- c("getVarCov.PREM", class(out))
  return(out)
}

#' @rdname getVarCov
#' @export
print.getVarCov.BPREM <- function(x, ...){

  cat("Random Effects Variance-Covariance Matrix\n")
  print(round(x$cov_mat,3), na.print="")
  cat("\n")
  cat("Random Effects Correlation Matrix\n")
  print(round(x$cor_mat,3), na.print="")
  cat("\n")

  invisible(x)
}

#' @rdname getVarCov
#' @export
print.getVarCov.CREM <- function(x, ...){

  cat("Random Effects Variance-Covariance Matrix for Individuals\n")
  print(round(x$ind_mat,3), na.print="")
  cat("\n")
  cat("Random Effects Variance-Covariance Matrix for Groups\n")
  print(round(x$grp_mat,3), na.print="")
  cat("\n")

  invisible(x)
}

#' @rdname getVarCov
#' @export
print.getVarCov.PREM <- function(x, ...){

  for(i in 1:length(x)){
    cat(paste0("Random Effects Variance-Covariance Matrix for Class ",i,"\n"))
    print(round(x[[i]],3), na.print="")
    cat("\n")
  }

  invisible(x)
}
