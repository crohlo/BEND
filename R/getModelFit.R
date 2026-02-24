#' Extract model fit
#'
#' @description
#' Extracts the model fit information from a fitted model of class "BPREM", "CREM", or "PREM".
#'
#' @param object An object of class "BPREM", "CREM", or "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a vector of the model fit information (deviance, pD, DIC).
#'
#' @author Corissa T. Rohloff
#'
#' @seealso [Bayes_PREM, Bayes_CREM, Bayes_PREM]
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get model fit
#' getModelFit(results_prem)
#'
#' @export
getModelFit <- function(object, ...) UseMethod(("getModelFit"))

#' @rdname getModelFit
#' @export
getModelFit.BPREM <- function(object, ...){

  out <- unlist(object$Model_Fit)
  class(out) <- c("getModelFit.BPREM", class(out))
  return(out)

}

#' @rdname getModelFit
#' @export
getModelFit.CREM <- function(object, ...){

  out <- unlist(object$Model_Fit)
  class(out) <- c("getModelFit.CREM", class(out))
  return(out)

}

#' @rdname getModelFit
#' @export
getModelFit.PREM <- function(object, ...){

  out <- unlist(object$Model_Fit)
  class(out) <- c("getModelFit.PREM", class(out))
  return(out)

}

#' @rdname getModelFit
#' @export
print.getModelFit <- function(object, ...){

  print(round(object,2), na.print="")
  cat("\n")

  invisible(object)
}

