#' Extract fixed effects parameter estimates
#'
#' @description
#' Extracts the fixed effects parameter estimates from a fitted model of class "BPREM", "CREM", or "PREM".
#'
#' @param object An object of class "BPREM", "CREM", or "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a vector or dataframe of fixed effects parameter values.
#'
#' @author Corissa T. Rohloff
#'
#' @seealso [Bayes_PREM, Bayes_CREM, Bayes_PREM]
#'
#' @examples
#' # load fitted model results
#' data(results_bprem)
#' # get fixed effects
#' getFixEf(results_bprem)
#'
#' @export
getFixEf <- function(object, ...) UseMethod(("getFixEf"))

#' @rdname getFixEf
#' @export
getFixEf.BPREM <- function(object, ...) {

  # determine number of parameters
  n_param <- 4
  # define parameters for labeling
  param_names <- c("Intercept", "Slope", "Change in Slope", "Changepoint")

  fix_eff_est <- matrix(object$Parameter_Estimates[1:(n_param*2),"Mean"], ncol=2)
  colnames(fix_eff_est) <- c(object$Call$y1_var, object$Call$y2_var)
  rownames(fix_eff_est) <- paste0(param_names)
  fix_eff_est <- as.data.frame(fix_eff_est)

  out <- list(fixef = fix_eff_est)
  class(out) <- c("getFixEf", class(out$fixef))
  return(out)
}

#' @rdname getFixEf
#' @export
getFixEf.CREM <- function(object, ...){

  # determine form
  form <- object$Functional_Form
  # determine number of parameters
  if(form=="linear")                           n_param <- 2
  if(form=="quadratic" | form=="exponential")  n_param <- 3
  if(form=="piecewise")                        n_param <- 4
  # define parameters for labeling
  if(form=="linear")       param_names <- c("Intercept", "Slope")
  if(form=="quadratic")    param_names <- c("Intercept", "Linear Slope", "Quadratic Slope")
  if(form=="exponential")  param_names <- c("Intercept", "Total Change", "Growth Rate")
  if(form=="piecewise")    param_names <- c("Intercept", "Slope", "Change in Slope", "Changepoint")

  fix_eff_est <- object$Parameter_Estimates[1:n_param,"Mean"]
  names(fix_eff_est) <- param_names

  out <- list(fixef = fix_eff_est)
  class(out) <- c("getFixEf", class(out$fixef))
  return(out)
}

#' @rdname getFixEf
#' @export
getFixEf.PREM <- function(object, ...){

  # determine number of classes
  n_class <- length(unique(object$Class_Information$class_membership))

  # determine number of changepoints in each class (based on final model results)
  changepoints <- c()
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    changepoints[i] <- which.max(object$Parameter_Estimates[[class_num]]$K_prob)-1
  }
  max_cp <- max(changepoints)

  # pull fixed effects for each class
  fixef_mat <- matrix(NA, nrow=2+2*max_cp, ncol=n_class)
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    cp_num <- paste0("K_", changepoints)[i]
    # fixed effects
    fixef_mat[1,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[1]
    fixef_mat[2,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[2]
    # if there are changepoints (beyond linear)
    if(changepoints[i]>0){
      for(k in 1:changepoints[i]){
        # fixed effects
        fixef_mat[2*k+1,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_mean[k]
        fixef_mat[2*k+2,i] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[k+2]
      }
    }
  }
  colnames(fixef_mat) <- paste0("Class ", 1:n_class)

  my_rownames <- rep(NA, 2+2*max_cp)
  my_rownames[1] = "Intercept"
  my_rownames[2] = "Slope"
  for(k in 1:max_cp){
    my_rownames[2*k+1] = paste0("Changepoint ", k)
    my_rownames[2*k+2] = paste0("Change in Slope ", k)
  }
  rownames(fixef_mat) <- my_rownames

  out <- list(fixef = fixef_mat)
  class(out) <- c("getFixEf", class(out$fixef))
  return(out)
}

#' @rdname getFixEf
#' @export
print.getFixEf <- function(object, ...){

  cat("Fixed Effects Parameters\n")
  print(object$fixef, na.print="")
  cat("\n")

  invisible(object)
}
