#' Simulated bivariate piecewise random effects model data
#'
#' A simulated dataset for the bivariate piecewise random effects model (BPREM).
#'
#' @format ## `SimData_BPREM`
#' A data frame with 210 rows and 4 columns:
#' \describe{
#'   \item{id}{Individual ids}
#'   \item{time}{Time}
#'   \item{y1}{Outcome 1}
#'   \item{y2}{Outcome 2}
#' }
"SimData_BPREM"

#' Simulated piecewise crossed random effects model data
#'
#' A simulated dataset for the piecewise crossed random effects model (PCREM).
#'
#' @format ## `SimData_PCREM`
#' A data frame with 210 rows and 4 columns:
#' \describe{
#'   \item{id}{Individual ids}
#'   \item{teacherid}{Teacher ids}
#'   \item{time}{Time}
#'   \item{y}{Outcome}
#' }
"SimData_PCREM"

#' Simulated piecewise random effects model data
#'
#' A simulated dataset for the piecewise random effects model (PREM).
#'
#' @format ## `SimData_PREM`
#' A data frame with 540 rows and 6 columns:
#' \describe{
#'   \item{id}{Individual ids}
#'   \item{time}{Time}
#'   \item{y}{Outcome}
#'   \item{class_pred_1}{First Class-Predictive Covariate}
#'   \item{class_pred_2}{Second Class-Predictive Covariate}
#'   \item{outcome_pred_1}{First Outcome-Predictive Covariate}
#' }
"SimData_PREM"

#' Bivariate piecewise random effects model results
#'
#' Results from fitting a bivariate piecewise random effects model (BPREM) to `SimData_BPREM` (returned by `Bayes_BPREM(...)`).
#'
#' @format ## `results_bprem`
#' A list object with elements:
#' \describe{
#'   \item{Call}{Function arguments}
#'   \item{Sample_Size}{Data set sample size}
#'   \item{Data}{Data used for model fitting}
#'   \item{Convergence}{Convergence status}
#'   \item{Model_Fit}{Model fit information}
#'   \item{Fitted_Values}{Fitted outcome values}
#'   \item{Parameter_Estimates}{Parameter estimates}
#'   \item{Random_Coefficients}{Random coefficient estimates}
#'   \item{Run_Time}{Model run time}
#' }
"results_bprem"

#' Piecewise crossed random effects model results
#'
#' Results from fitting a piecewise crossed random effects model (PCREM) to `SimData_PCREM` (returned by `Bayes_CREM(...)`).
#'
#' @format ## `results_pcrem`
#' A list object with elements:
#' \describe{
#'   \item{Call}{Function arguments}
#'   \item{Sample_Size}{Data set sample size}
#'   \item{Data}{Data used for model fitting}
#'   \item{Convergence}{Convergence status}
#'   \item{Model_Fit}{Model fit information}
#'   \item{Fitted_Values}{Fitted outcome values}
#'   \item{Parameter_Estimates}{Parameter estimates}
#'   \item{Random_Coefficients}{Random coefficient estimates}
#'   \item{Run_Time}{Model run time}
#' }
"results_pcrem"

#' Piecewise random effects model results
#'
#' Results from fitting a piecewise random effects model (PREM) to `SimData_PREM` (returned by `Bayes_PREM(...)`).
#'
#' @format ## `results_prem`
#' A list object with elements:
#' \describe{
#'   \item{Call}{Function arguments}
#'   \item{Sample_Size}{Data set sample size}
#'   \item{Data}{Data used for model fitting}
#'   \item{Convergence}{Convergence status}
#'   \item{Model_Fit}{Model fit information}
#'   \item{Fitted_Values}{Fitted outcome values}
#'   \item{Parameter_Estimates}{Parameter estimates}
#'   \item{Random_Coefficients}{Random coefficient estimates}
#'   \item{Run_Time}{Model run time}
#' }
"results_prem"
