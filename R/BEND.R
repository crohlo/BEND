#' Simulated data for a PREM + Extensions
#'
#' Simulated data for a piecewise random effects model (PREM) and useful extensions (CI-PREM, PREMM, CI-PREMM) with 18 timepoints collected on 30 individuals.
#'
#' \itemize{
#'   \item `id` ID for each individual.
#'   \item `time` Timepoints for each individual.
#'   \item `y` Outcome.
#'   \item `class_pred_1` First class predictive covariate (time-invariant).
#'   \item `class_pred_2` Second class predictive covariate (time-invariant).
#'   \item `outcome_pred_1` Outcome predictive covariate (time-varying).
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SimData_PREM
#' @usage data(SimData_PREM)
#' @format A data frame with 540 rows and 6 variables.
NULL

#' Simulated data for a PCREM
#'
#' Simulated data for a piecewise crossed random effects model (PCREM) with 7 timepoints collected on 30 individuals.
#'
#' \itemize{
#'   \item `id` ID for each individual.
#'   \item `teacherid` ID for each teacher.
#'   \item `time` Timepoints for each individual.
#'   \item `y` Outcome.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SimData_PCREM
#' @usage data(SimData_PCREM)
#' @format A data frame with 210 rows and 4 variables.
NULL

#' Simulated data for a BPREM
#'
#' Simulated data for a bivariate piecewise random effects model (BPREM) with 7 timepoints collected on 30 individuals.
#'
#' \itemize{
#'   \item `id` ID for each individual.
#'   \item `time` Timepoints for each individual.
#'   \item `y1` Outcome 1.
#'   \item `y2` Outcome 2.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name SimData_BPREM
#' @usage data(SimData_BPREM)
#' @format A data frame with 210 rows and 4 variables.
NULL

#' Fitted results for a PREM
#'
#' Fitted results for a piecewise random effects model (PREM) using `SimData_PREM`. Included to demonstrate the use of `plot_BEND()` and `summary.PREM()`.
#'
#' @docType data
#' @keywords datasets
#' @name results_prem
#' @usage data(results_prem)
#' @format A list (an object of class `PREM`) with fitted model results.
NULL

#' Fitted results for a PCREM
#'
#' Simulated data for a piecewise crossed random effects model (PCREM) using `SimData_CREM`. Included to demonstrate the use of `summary.CREM()`.
#'
#' @docType data
#' @keywords datasets
#' @name results_pcrem
#' @usage data(results_pcrem)
#' @format A list (an object of class `CREM`) with fitted model results.
NULL

#' Fitted results for a BPREM
#'
#' Simulated data for a bivariate piecewise random effects model (BPREM) using `SimData_BPREM`. Included to demonstrate the use of `summary.BPREM()`.
#'
#' @docType data
#' @keywords datasets
#' @name results_bprem
#' @usage data(results_bprem)
#' @format A list (an object of class `BPREM`) with fitted model results.
NULL
