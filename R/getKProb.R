#' Extract changepoint probabilities
#'
#' @description
#' Extracts the K (number of changepoints) probabilities from a fitted model of class "PREM".
#'
#' @param x An object of class "PREM".
#' @param ... Additional arguments.
#'
#' @returns Returns a list of the posterior probabilities for each possible number of changepoints, for each class.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # get changepoint probabilities
#' getKProb(results_prem)
#'
#' @export
getKProb <- function(x, ...) UseMethod(("getKProb"))

#' @rdname getKProb
#' @export
getKProb.PREM <- function(x, ...){

  # determine number of classes
  n_class <- length(unique(x$Class_Information$class_membership))

  k_probs <- vector("list", n_class)
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    k_probs[[i]] <- x$Parameter_Estimates[[class_num]]$K_prob
  }
  names(k_probs) <- paste0("Class ", 1:n_class)

  out <- k_probs
  class(out) <- c("getKProb.PREM", class(out))
  return(out)

}

#' @rdname getKProb
#' @export
print.getKProb.PREM <- function(x, ...){

  cat("K (Number of Changepoints) Probabilities\n\n")
  for(i in 1:length(x)){
    cat(paste0("Class ",i,"\n"))
    print(x[[i]], na.print="")
    cat("\n")
  }

  invisible(x)
}

