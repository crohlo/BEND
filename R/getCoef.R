#' Extract random coefficients
#'
#' @description
#' Extracts the random coefficients from a fitted model of class "BPREM", "CREM", or "PREM".
#'
#' @param x An object of class "BPREM", "CREM", or "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a data frame of the random coefficients for each individual/group.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get random coefficients
#' getCoef(results_prem)
#'
#' @export
getCoef <- function(x, ...) UseMethod(("getCoef"))

#' @rdname getCoef
#' @export
getCoef.BEND <- function(x, ...){

  out <- x$Random_Coefficients$rancoef
  class(out) <- c("getCoef.BEND", class(out))
  return(out)

}

#' @rdname getCoef
#' @export
print.getCoef <- function(x, ...){

  print(x, na.print="")
  cat("\n")

  invisible(x)
}

