#' Extract model fit
#'
#' @description
#' Extracts the model fit information from a fitted model of class "BPREM", "CREM", or "PREM".
#'
#' @param x An object of class "BPREM", "CREM", or "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a vector of the model fit information (deviance, pD, DIC).
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get model fit
#' getModelFit(results_prem)
#'
#' @export
getModelFit <- function(x, ...) UseMethod(("getModelFit"))

#' @rdname getModelFit
#' @export
getModelFit.BEND <- function(x, ...){

  out <- list(mod_fit = as.data.frame(x$Model_Fit))
  class(out) <- c("getModelFit.BEND", class(out))
  return(out)

}

#' @rdname getModelFit
#' @export
print.getModelFit.BEND <- function(x, ...){

  print(round(x$mod_fit,2), na.print="", row.names = FALSE)
  cat("\n")

  invisible(x)
}

