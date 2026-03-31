#' Extract class probabilities
#'
#' @description
#' Extracts the posterior class probabilities for each individual from a fitted model of class "PREM".
#'
#' @param x An object of class "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a list of the posterior class probabilities for each individual.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get class probabilities
#' getClassProb(results_prem)
#'
#' @export
getClassProb <- function(x, ...) UseMethod(("getClassProb"))

#' @rdname getClassProb
#' @export
getClassProb.PREM <- function(x, ...){

  # determine number of classes
  n_class <- length(unique(x$Class_Information$class_membership))

  # detect covariates
  n_cov_op <- length(x$Parameter_Estimates$outcome_predictive_covariates)
  n_cov_cp <- length(x$Parameter_Estimates$class_predictive_covariates)

  if(n_class==1){
    individ_class_info <- cbind(names(x$Class_Information$class_membership),
                                x$Class_Information$class_membership,
                                x$Class_Information$individ_class_probability)
    individ_class_info <- as.data.frame(individ_class_info)
    names(individ_class_info) <- c(x$Call$id_var, "class_membership", "probability_class_1")
    class_prob <- x$Class_Information$unconditional_class_probability
  }
  if(n_class>1){
    if(n_cov_cp==0){
      individ_class_info <- cbind(names(x$Class_Information$class_membership),
                                  x$Class_Information$class_membership,
                                  x$Class_Information$individ_class_probability)
      individ_class_info <- as.data.frame(individ_class_info)
      names(individ_class_info) <- c(x$Call$id_var, "class_membership", paste0("probability_class_",1:n_class))
      class_prob <- x$Class_Information$unconditional_class_probability
    }
    if(n_cov_cp>0){
      individ_class_info <- cbind(names(x$Class_Information$class_membership),
                                  x$Class_Information$class_membership,
                                  x$Class_Information$individ_class_probability,
                                  x$Class_Information$unconditional_class_probability)
      individ_class_info <- as.data.frame(individ_class_info)
      names(individ_class_info) <- c(x$Call$id_var, "class_membership", paste0("probability_class_",1:n_class), "probability_class_ref")
      class_prob <- NULL
    }
  }

  out <- list(class_prob = class_prob,
              individ_class_info = individ_class_info)
  class(out) <- c("getClassProb.PREM", class(out))
  return(out)

}

#' @rdname getClassProb
#' @export
print.getClassProb.PREM <- function(x, ...){

  if(!is.null(x$class_prob)){
    cat("Population Class Probabilities\n")
    print(round(x$class_prob,2), na.print="")
    cat("\n")
  }

  cat("Individual Class Information\n")
  print(x$individ_class_info, na.print="")
  cat("\n")

  invisible(x)
}

