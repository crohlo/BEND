#' Extract random coefficients
#'
#' @description
#' Extracts the random coefficients from a fitted model of class "BPREM", "CREM", or "PREM".
#'
#' @param object An object of class "BPREM", "CREM", or "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a data frame of the random coefficients for each individual/group.
#'
#' @author Corissa T. Rohloff
#'
#' @seealso [Bayes_PREM, Bayes_CREM, Bayes_PREM]
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get random coefficients
#' getCoef(results_prem)
#'
#' @export
getCoef <- function(object, ...) UseMethod(("getCoef"))

#' @rdname getCoef
#' @export
getCoef.BPREM <- function(object, ...){

  out <- object$Random_Coefficients$rancoef
  class(out) <- c("getCoef.BPREM", class(out))
  return(out)

}

#' @rdname getCoef
#' @export
getCoef.CREM <- function(object, ...){

  out <- object$Random_Coefficients$rancoef
  class(out) <- c("getCoef.CREM", class(out))
  return(out)

}

#' @rdname getCoef
#' @export
getCoef.PREM <- function(object, ...){

  out <- object$Random_Coefficients$rancoef
  class(out) <- c("getCoef.PREM", class(out))
  return(out)

}

#' @rdname getCoef
#' @export
print.getCoef <- function(object, ...){

  print(object, na.print="")
  cat("\n")

  invisible(object)
}

