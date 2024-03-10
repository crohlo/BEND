#' Simulated data for a PREM
#'
#' Simulated data for a piecewise random effects model (PREM) with 10 timepoints collected on 20 individuals. The variables are as follows:
#'
#' \itemize{
#'   \item id - id for each individual
#'   \item time - timepoints for each individual
#'   \item y - outcome
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SimData_PREM
#' @usage data(SimData_PREM)
#' @format A data frame with 200 rows and 3 variables
NULL

#' Simulated data for a CI-PREMM
#'
#' Simulated data for a covariate-influenced piecewise random effects mixture model (CI-PREMM) with 18 timepoints collected on 30 individuals. The variables are as follows:
#'
#' \itemize{
#'   \item id - id for each individual
#'   \item time - timepoints for each individual
#'   \item y - outcome
#'   \item class_pred_1 - first class predictive covariate (time-invariant)
#'   \item class_pred_2 - second class predictive covariate (time-invariant)
#'   \item outcome_pred_1 - outcome predictive covariate (time-varying)
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SimData_CIPREMM
#' @usage data(SimData_CIPREMM)
#' @format A data frame with 540 rows and 6 variables
NULL

#' Simulated data for a PCREM
#'
#' Simulated data for a piecewise crossed random effects model (PCREM) with 7 timepoints collected on 30 individuals. The variables are as follows:
#'
#' \itemize{
#'   \item id - id for each individual
#'   \item teacherid - id for each teacher
#'   \item time - timepoints for each individual
#'   \item y - outcome
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SimData_PCREM
#' @usage data(SimData_PCREM)
#' @format A data frame with 210 rows and 4 variables
NULL

#' Simulated data for a BPREM
#'
#' Simulated data for a bivariate piecewise random effects model (BPREM) with 7 timepoints collected on 30 individuals. The variables are as follows:
#'
#' \itemize{
#'   \item id - id for each individual
#'   \item time - timepoints for each individual
#'   \item y1 - outcome 1
#'   \item y2 - outcome 2
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SimData_BPREM
#' @usage data(SimData_BPREM)
#' @format A data frame with 210 rows and 4 variables
NULL
