#' Extract fitted values
#'
#' @description
#' Extracts the individual fitted values from a fitted model of class "BPREM", "CREM", or "PREM".
#'
#' @param object An object of class "BPREM", "CREM", or "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a data frame of the fitted values for each individual and timepoint.
#'
#' @author Corissa T. Rohloff
#'
#' @seealso [Bayes_PREM, Bayes_CREM, Bayes_PREM]
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get fitted values
#' getFitted(results_prem)
#'
#' @export
getFitted <- function(object, ...) UseMethod(("getFitted"))

#' @rdname getFitted
#' @export
getFitted.BPREM <- function(object, ...){

  out <- object$Fitted_Values
  class(out) <- c("getFitted.BPREM", class(out))
  return(out)

}

#' @rdname getFitted
#' @export
getFitted.CREM <- function(object, ...){

  out <- object$Fitted_Values
  class(out) <- c("getFitted.CREM", class(out))
  return(out)

}

#' @rdname getFitted
#' @export
getFitted.PREM <- function(object, ...){

  out <- object$Fitted_Values
  class(out) <- c("getFitted.PREM", class(out))
  return(out)

}

#' @rdname getFitted
#' @export
print.getFitted <- function(object, ...){

  print(object, na.print="")
  cat("\n")

  invisible(object)
}

