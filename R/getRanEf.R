#' Extract random effects
#'
#' @description
#' Extracts the random effects from a fitted model of class "BPREM" or "CREM".
#'
#' @param object An object of class "BPREM" or "CREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a list of the random effects for each individual/group.
#'
#' @author Corissa T. Rohloff
#'
#' @seealso [Bayes_BPREM, Bayes_CREM]
#'
#' @examples
#' # load fitted model results
#' data(results_bprem)
#' # get random effects
#' getRanEf(results_bprem)
#'
#' @export
getRanEf <- function(object, ...) UseMethod(("getRanEf"))

#' @rdname getRanEf
#' @export
getRanEf.BPREM <- function(object, ...){

  out <- list(ranef = object$Random_Coefficients$ranef_b)
  class(out) <- c("getRanEf.BPREM", class(out))
  return(out)

}

#' @rdname getRanEf
#' @export
getRanEf.CREM <- function(object, ...){

  out <- list(ranef_b = object$Random_Coefficients$ranef_b,
              ranef_g = object$Random_Coefficients$ranef_g)
  class(out) <- c("getRanEf.CREM", class(out))
  return(out)

}

#' @rdname getRanEf
#' @export
print.getRanEf <- function(object, ...){

  print(object, na.print="")
  cat("\n")

  invisible(object)
}

