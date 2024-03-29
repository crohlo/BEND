% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bayes_BPREM.R
\name{Bayes_BPREM}
\alias{Bayes_BPREM}
\title{Bayesian Bivariate Piecewise Random Effects Model (BPREM)}
\usage{
Bayes_BPREM(
  data,
  id_var,
  time_var,
  y1_var,
  y2_var,
  iters_adapt = 5000,
  iters_burn_in = 1e+05,
  iters_sampling = 50000,
  thin = 15,
  save_full_chains = FALSE,
  save_conv_chains = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{data}{Data frame in long format, where each row describes a measurement occasion for a given individual. It is assumed that each individual has the same number of assigned timepoints (a.k.a., rows). There can be missingness in the outcomes (\code{y1_var} and \code{ y2_var}), but there cannot be missingness in time (\code{time_var}).}

\item{id_var}{Name of column that contains ids for individuals with repeated measures in a longitudinal dataset.}

\item{time_var}{Name of column that contains the time variable. This column cannot contain any missing values.}

\item{y1_var}{Name of column that contains the first outcome variable. Missing values should be denoted by NA.}

\item{y2_var}{Name of column that contains the second outcome variable. Missing values should be denoted by NA.}

\item{iters_adapt}{(optional) Number of iterations for adaptation of jags model (default = 5000).}

\item{iters_burn_in}{(optional) Number of iterations for burn-in (default = 100000).}

\item{iters_sampling}{(optional) Number of iterations for posterior sampling (default = 50000).}

\item{thin}{(optional) Thinning interval for posterior sampling (default = 15).}

\item{save_full_chains}{Logical indicating whether the MCMC chains from rjags should be saved (default = FALSE). Note, this should not be used regularly as it will result in an object with a large file size.}

\item{save_conv_chains}{Logical indicating whether the MCMC chains from rjags should be saved but only for the parameters monitored for convergence (default = FALSE). This would be useful for plotting traceplots for relevant model parameters to evaluate convergence behavior. Note, this should not be used regularly as it will result in an object with a large file size.}

\item{verbose}{Logical controlling whether progress messages/bars are generated (default = TRUE).}
}
\value{
A list (an object of class \code{BPREM}) with elements:
\item{Convergence}{Potential scale reduction factor (PSRF) for each parameter (\code{parameter_psrf}), Gelman multivariate scale reduction factor (\code{multivariate_psrf}), and mean PSRF (\code{mean_psrf}) to assess model convergence.}
\item{Model_Fit}{Deviance (\code{deviance}), effective number of parameters (\code{pD}), and Deviance information criterion (\code{dic}) to assess model fit.}
\item{Fitted_Values}{Vector giving the fitted value at each timepoint for each individual (same length as long data).}
\item{Parameter_Estimates}{Data frame with posterior mean and 95\% credible intervals for each model parameter.}
\item{Run_Time}{Total run time for model fitting.}
\item{Full_MCMC_Chains}{If save_full_chains=TRUE, raw MCMC chains from rjags.}
\item{Convergence_MCMC_Chains}{If save_conv_chains=TRUE, raw MCMC chains from rjags but only for the parameters monitored for convergence.}
}
\description{
Estimates a Bayesian bivariate piecewise random effects models (BPREM) for longitudinal data with two interrelated outcomes. See Peralta et al. (2022) for more details.
}
\details{
For more information on the model equation and priors implemented in this function, see Peralta et al. (2022).
}
\examples{
\donttest{
# load simulated data
data(SimData_BPREM)
# plot observed data
plot_BEND(data = SimData_BPREM,
          id_var = "id",
          time_var = "time",
          y_var = "y1",
          y2_var = "y2")
# fit Bayes_BPREM()
results_bprem <- Bayes_BPREM(data = SimData_BPREM,
                             id_var = "id",
                             time_var = "time",
                             y1_var = "y1",
                             y2_var = "y2")
# result summary
summary(results_bprem)
# plot fitted results
plot_BEND(data = SimData_BPREM,
          id_var = "id",
          time_var = "time",
          y_var = "y1",
          y2_var = "y2",
          results = results_bprem)
}

}
\references{
Peralta, Y., Kohli, N., Lock, E. F., & Davison, M. L. (2022). Bayesian modeling of associations in bivariate piecewise linear mixed-effects models. Psychological Methods, 27(1), 44–64. https://doi.org/10.1037/met0000358
}
\author{
Corissa T. Rohloff, Yadira Peralta
}
