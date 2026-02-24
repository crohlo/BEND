#' Extract class probabilities
#'
#' @description
#' Extracts the class probabilities from a fitted model of class "PREM".
#'
#' @param object An object of class "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a list of the class probabilities.
#'
#' @author Corissa T. Rohloff
#'
#' @seealso [Bayes_PREM]
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get class probabilities
#' getClassProb(results_prem)
#'
#' @export
getClassProb <- function(object, ...) UseMethod(("getClassProb"))

#' @rdname getClassProb
#' @export
getClassProb.PREM <- function(object, ...){

  # determine number of classes
  n_class <- length(unique(object$Class_Information$class_membership))

  # detect covariates
  n_cov_op <- length(object$Parameter_Estimates$outcome_predictive_covariates)
  n_cov_cp <- length(object$Parameter_Estimates$class_predictive_covariates)

  if(n_class==1){
    individ_class_info <- cbind(names(object$Class_Information$class_membership),
                                object$Class_Information$class_membership,
                                object$Class_Information$individ_class_probability)
    individ_class_info <- as.data.frame(individ_class_info)
    names(individ_class_info) <- c(object$Call$id_var, "class_membership", "probability_class_1")
    class_prob <- object$Class_Information$unconditional_class_probability
  }
  if(n_class>1){
    if(n_cov_cp==0){
      individ_class_info <- cbind(names(object$Class_Information$class_membership),
                                  object$Class_Information$class_membership,
                                  object$Class_Information$individ_class_probability)
      individ_class_info <- as.data.frame(individ_class_info)
      names(individ_class_info) <- c(object$Call$id_var, "class_membership", paste0("probability_class_",1:n_class))
      class_prob <- object$Class_Information$unconditional_class_probability
    }
    if(n_cov_cp>0){
      individ_class_info <- cbind(names(object$Class_Information$class_membership),
                                  object$Class_Information$class_membership,
                                  object$Class_Information$individ_class_probability,
                                  object$Class_Information$unconditional_class_probability)
      individ_class_info <- as.data.frame(individ_class_info)
      names(individ_class_info) <- c(object$Call$id_var, "class_membership", paste0("probability_class_",1:n_class), "probability_class_ref")
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
print.getClassProb.PREM <- function(object, ...){

  if(!is.null(object$class_prob)){
    cat("Population Class Probabilities\n")
    print(round(object$class_prob,2), na.print="")
    cat("\n")
  }

  cat("Individual Class Information\n")
  print(object$individ_class_info, na.print="")
  cat("\n")

  invisible(object)
}

