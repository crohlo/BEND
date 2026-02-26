#' Extract random effects
#'
#' @description
#' Extracts the random effects from a fitted model of class "BPREM" or "CREM".
#'
#' @param x An object of class "BPREM" or "CREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a list of the random effects for each individual/group.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_bprem)
#' # get random effects
#' getRanEf(results_bprem)
#'
#' @export
getRanEf <- function(x, ...) UseMethod(("getRanEf"))

#' @rdname getRanEf
#' @export
getRanEf.BPREM <- function(x, ...){

  out <- list(ranef = x$Random_Coefficients$ranef_b)
  class(out) <- c("getRanEf.BPREM", class(out))
  return(out)

}

#' @rdname getRanEf
#' @export
getRanEf.CREM <- function(x, ...){

  out <- list(ranef_b = x$Random_Coefficients$ranef_b,
              ranef_g = x$Random_Coefficients$ranef_g)
  class(out) <- c("getRanEf.CREM", class(out))
  return(out)

}

#' @rdname getRanEf
#' @export
print.getRanEf <- function(x, ...){

  print(x, na.print="")
  cat("\n")

  invisible(x)
}

