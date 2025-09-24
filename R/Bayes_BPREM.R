#' Bayesian Bivariate Piecewise Random Effects Model (BPREM)
#'
#' @description Estimates a Bayesian bivariate piecewise random effects models (BPREM) for longitudinal data with two interrelated outcomes. See Peralta et al. (2022) for more details.
#'
#' @param data Data frame in long format, where each row describes a measurement occasion for a given individual. It is assumed that each individual has the same number of assigned timepoints (a.k.a., rows). There can be missingness in the outcomes (`y1_var` and ` y2_var`), but there cannot be missingness in time (`time_var`).
#' @param id_var Name of column that contains ids for individuals with repeated measures in a longitudinal dataset.
#' @param time_var Name of column that contains the time variable. This column cannot contain any missing values.
#' @param y1_var Name of column that contains the first outcome variable. Missing values should be denoted by NA.
#' @param y2_var Name of column that contains the second outcome variable. Missing values should be denoted by NA.
#' @param iters_adapt (optional) Number of iterations for adaptation of jags model (default = 5000).
#' @param iters_burn_in (optional) Number of iterations for burn-in (default = 100000).
#' @param iters_sampling (optional) Number of iterations for posterior sampling (default = 50000).
#' @param thin (optional) Thinning interval for posterior sampling (default = 15).
#' @param save_full_chains Logical indicating whether the MCMC chains from rjags should be saved (default = FALSE). Note, this should not be used regularly as it will result in an object with a large file size.
#' @param save_conv_chains Logical indicating whether the MCMC chains from rjags should be saved but only for the parameters monitored for convergence (default = FALSE). This would be useful for plotting traceplots for relevant model parameters to evaluate convergence behavior. Note, this should not be used regularly as it will result in an object with a large file size.
#' @param verbose Logical controlling whether progress messages/bars are generated (default = TRUE).
#'
#' @returns A list (an object of class `BPREM`) with elements:
#' \item{Convergence}{Potential scale reduction factor (PSRF) for each parameter (`parameter_psrf`), Gelman multivariate scale reduction factor (`multivariate_psrf`), and mean PSRF (`mean_psrf`) to assess model convergence.}
#' \item{Model_Fit}{Deviance (`deviance`), effective number of parameters (`pD`), and Deviance information criterion (`dic`) to assess model fit.}
#' \item{Fitted_Values}{Vector giving the fitted value at each timepoint for each individual (same length as long data).}
#' \item{Parameter_Estimates}{Data frame with posterior mean and 95% credible intervals for each model parameter.}
#' \item{Run_Time}{Total run time for model fitting.}
#' \item{Full_MCMC_Chains}{If save_full_chains=TRUE, raw MCMC chains from rjags.}
#' \item{Convergence_MCMC_Chains}{If save_conv_chains=TRUE, raw MCMC chains from rjags but only for the parameters monitored for convergence.}
#'
#' @details
#' For more information on the model equation and priors implemented in this function, see Peralta et al. (2022).
#'
#' @author Corissa T. Rohloff, Yadira Peralta
#'
#' @references Peralta, Y., Kohli, N., Lock, E. F., & Davison, M. L. (2022). Bayesian modeling of associations in bivariate piecewise linear mixed-effects models. Psychological Methods, 27(1), 44–64. https://doi.org/10.1037/met0000358
#'
#' @examples
#' \donttest{
#' # load simulated data
#' data(SimData_BPREM)
#' # plot observed data
#' plot_BEND(data = SimData_BPREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y1",
#'           y2_var = "y2")
#' # fit Bayes_BPREM()
#' results_bprem <- Bayes_BPREM(data = SimData_BPREM,
#'                              id_var = "id",
#'                              time_var = "time",
#'                              y1_var = "y1",
#'                              y2_var = "y2")
#' # result summary
#' summary(results_bprem)
#' # plot fitted results
#' plot_BEND(data = SimData_BPREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y1",
#'           y2_var = "y2",
#'           results = results_bprem)
#' }
#'
#' @import stats
#'
#' @export
Bayes_BPREM <- function(data,
                        id_var, time_var, y1_var, y2_var,
                        iters_adapt=5000, iters_burn_in=100000, iters_sampling=50000, thin=15,
                        save_full_chains=FALSE, save_conv_chains=FALSE,
                        verbose=TRUE){

  # START OF SETUP ----

  ## Initial data check
  data <- as.data.frame(data)

  ## Control progress messages/bars
  if(verbose) progress_bar = "text"
  if(!verbose) progress_bar = "none"

  ## Start tracking run time
  run_time_total_start <- Sys.time()

  ## Load module to compute DIC
  suppressMessages(rjags::load.module('dic'))

  ## Set number of chains - will use 3 chains across all models
  n_chains <- 3

  # Reshape data for modeling ----

  ## outcome data - matrix form
  y1 <- reshape(data[,c(id_var, time_var, y1_var)],
                idvar=id_var,
                timevar=time_var,
                direction='wide')
  y1 <- unname(as.matrix(y1[,names(y1)!=id_var]))

  y2 <- reshape(data[,c(id_var, time_var, y2_var)],
                idvar=id_var,
                timevar=time_var,
                direction='wide')
  y2 <- unname(as.matrix(y2[,names(y2)!=id_var]))

  ## JAGS model assumes the two outcome variables are in an array
  y <- array(c(y1, y2), dim = c(nrow(y1),ncol(y1),2))

  ## time data - matrix form
  ## should be the same dimensions as y1 and y2
  x <- matrix(data[,c(time_var)],
              byrow=TRUE,
              nrow=dim(y1)[1],
              ncol=dim(y1)[2])

  # Define relevant variables ----
  n_subj <- dim(y)[1]
  n_time <- dim(y)[2]

  max_time <- max(x)
  min_time <- min(x)
  min_cp_mean <- unique(sort(x))[2]
  max_cp_mean <- unique(sort(x, decreasing=TRUE))[2]

  mean_cp <- (max_time-min_time)/2          # mean of changepoint
  prec_cp <- 1/((max_time-min_time)/4)^2    # precision of changepoint
  bound_e1 <- min(apply(y[,,1], 2, var, na.rm = TRUE))  # bound of error variance 1
  bound_e2 <- min(apply(y[,,2], 2, var, na.rm = TRUE))  # bound of error variance 2
  beta_mean_zero <- rep(0,8)  # Vector of zeros to center the distribution of raw random-effects
  omega_b <- diag(8)     # Scale matrix for the Wishart distribution

  # Error messages for input ----

  ## data format warnings
  if(sum(is.na(x)>0)) stop('Columns for t_var cannot have NA values (but y_var can)')
  if(min(x)!=0) warning('Prior assumes first time point measured at x=0')

  # FULL MODEL ----

  ## Specify full JAGS model
  full_spec <- textConnection(bivariate_pw)

  ## Variables to extract from full model
  param_recovery_full <- c('beta', 'beta_mean',
                           'var_b', 'cov_b',
                           'cor_b', 'rho',
                           'mu_y', 'sigma2_error',
                           'for_conv', 'pD', 'deviance')

  ## Create data list for JAGS model
  data_list <- list('y' = y, 'x' = x,
                    'n_subj' = n_subj, 'n_time' = n_time,
                    'mean_cp' = mean_cp, 'prec_cp' = prec_cp,
                    'min_time' = min_time, 'max_time' = max_time,
                    'bound_e1' = bound_e1, 'bound_e2' = bound_e2,
                    'beta_mean_zero' = beta_mean_zero,
                    'omega_b' = omega_b)

  ## Set Seeds for Reproducible Code
  initial_vals <- vector('list', n_chains)
  seeds <- sample(1:10, n_chains, replace = FALSE) # randomly generate non-repeating integers
  # sensitive to set.seed()
  for(i in 1:n_chains){
    initial_vals[[i]]$.RNG.name <- "base::Wichmann-Hill" # arbitrary
    initial_vals[[i]]$.RNG.seed <- seeds[i]
  }

  ## Create JAGS model
  if(verbose) cat('Calibrating MCMC...\n')
  full_model <- rjags::jags.model(full_spec,
                                  data = data_list,
                                  inits = initial_vals,
                                  n.chains = n_chains,
                                  n.adapt = iters_adapt,
                                  quiet = TRUE)

  # burn-in
  if(verbose) cat('Running burn-in...\n')
  update(full_model, iters_burn_in, progress.bar=progress_bar)

  # sampling
  if(verbose) cat('Collecting samples...\n')
  full_out <- rjags::jags.samples(full_model,
                                  variable.names=param_recovery_full,
                                  n.iter=iters_sampling,
                                  thin=thin,
                                  progress.bar=progress_bar)

  # Compiling Full Results -----

  ## Convergence

  # Setup the parameter vector
  beta_mean_param <- paste0('beta_', rep(1:2,e=4), c(0:2, 'cp'), '_mean')
  beta_var_param <- paste0('var_b_', rep(1:2,e=4), c(0:2, 'cp'))
  base_cov <- c('11_10',
                '12_10',  '12_11',
                '1cp_10', '1cp_11', '1cp_12',
                '20_10',  '20_11',  '20_12',  '20_1cp',
                '21_10',  '21_11',  '21_12',  '21_1cp',  '21_20',
                '22_10',  '22_11',  '22_12',  '22_1cp',  '22_20',  '22_21',
                '2cp_10', '2cp_11', '2cp_12', '2cp_1cp', '2cp_20', '2cp_21', '2cp_22')
  beta_cov_param <- paste0("cov_b_", base_cov)
  beta_corr_param <- paste0("corr_b_", base_cov)
  error_param <- c("error_var_11", "error_var_22", "error_cov_12")
  param_names <- c(beta_mean_param, beta_var_param, beta_cov_param, error_param,
                   beta_corr_param, "error_corr_12")

  mcmc_list <- coda::as.mcmc.list(full_out$for_conv)
  for(i in 1:3){colnames(mcmc_list[[i]]) <- param_names}

  gelman_msrf <- coda::gelman.diag(mcmc_list[,1:47]) # individual parameter psrf, and multivariate psrf
  # beta_corr_param and error_corr_12 are redundant parameters, thus they were removed when assessing convergence
  parameter_psrf <- data.frame(point_est = gelman_msrf$psrf[,1],
                               upper_ci = gelman_msrf$psrf[,2])
  multivariate_psrf <- gelman_msrf$mpsrf
  mean_psrf <- mean(parameter_psrf$point_est) # mean psrf across parameters
  convergence <- list(parameter_psrf=parameter_psrf, multivariate_psrf=multivariate_psrf, mean_psrf=mean_psrf)

  ## Model Fit
  deviance <- mean(full_out$deviance)
  pD <- mean(full_out$pD)
  dic <- deviance + pD
  model_fit <- list(deviance=deviance, pD=pD, dic=dic)

  ## Fitted Values
  y_mean <- summary(full_out$mu_y, FUN='mean')[[1]]

  ## Parameter Estimates
  sum_mcmc <- summary(mcmc_list) # parameter estimates
  ## Parameter Estimates
  sum_mcmc <- summary(mcmc_list) # parameter estimates
  param_est <- data.frame(Mean = sum_mcmc$statistics[,1],
                          CI_Lower = sum_mcmc$quantiles[,1],
                          CI_Upper = sum_mcmc$quantiles[,5])

  ## Stop tracking run time
  run_time_total_end <- Sys.time()
  run_time_total <- run_time_total_end - run_time_total_start

  my_results <- list('Convergence' = convergence,
                     'Model_Fit' = model_fit,
                     'Fitted_Values'=y_mean,
                     'Parameter_Estimates'=param_est,
                     'Run_Time'=format(run_time_total))
  if(save_full_chains==TRUE){my_results$Full_MCMC_Chains=full_out}
  if(save_conv_chains==TRUE){my_results$Convergence_MCMC_Chains=mcmc_list[,1:47]}
  class(my_results) <- 'BPREM'
  return(my_results)
}

# Model ----

## Bivariate Piecewise ----

bivariate_pw <- "model{
  for(i in 1:n_subj){
    for(j in 1:n_time){
      y[i,j,1:2] ~ dmnorm(mu_y[i,j,1:2], tau_error[1:2,1:2])
      mu_y[i,j,1] <- beta[i,1] + beta[i,2]*x[i,j] + beta[i,3]*(max(0, x[i,j]-beta[i,4]))
      mu_y[i,j,2] <- beta[i,5] + beta[i,6]*x[i,j] + beta[i,7]*(max(0, x[i,j]-beta[i,8]))
    }
  }

  ## Level 1 precision and error variance
  tau_error[1:2,1:2] <- inverse(sigma2_e[,])
  sigma2_e[1,1] <- sigma2_e1
  sigma2_e[2,2] <- sigma2_e2
  sigma2_e[2,1] <- rho*sqrt(sigma2_e1)*sqrt(sigma2_e2)
  sigma2_e[1,2] <- sigma2_e[2,1]
  sigma2_e1 ~ dunif(0,bound_e1)
  sigma2_e2 ~ dunif(0,bound_e2)
  rho ~ dunif(-1,1)

  ## Distribution of random-effects and error variance
  for (i in 1:n_subj){
    for(k in 1:8){
      beta[i,k] <- beta_mean[k] + beta_rand[i,k]
      beta_rand[i,k] <- scale_c[k]*beta_raw[i,k]
    }
    beta_raw[i,1:8] ~ dmnorm(beta_mean_zero[1:8], tau_b_raw[1:8,1:8])
  }

  ## Priors for fixed-effects
  beta_mean[1] ~ dnorm(0, 0.0001)
  beta_mean[2] ~ dnorm(0, 0.0001)
  beta_mean[3] ~ dnorm(0, 0.0001)
  beta_mean[4] ~ dnorm(mean_cp, prec_cp)T(min_time,max_time)
  beta_mean[5] ~ dnorm(0, 0.0001)
  beta_mean[6] ~ dnorm(0, 0.0001)
  beta_mean[7] ~ dnorm(0, 0.0001)
  beta_mean[8] ~ dnorm(mean_cp, prec_cp)T(min_time,max_time)

  ## Priors for scaling constants
  for(k in 1:8){
    scale_c[k] ~ dunif(0, 100)
  }

  ## Prior for covariance matrix of random-effects
  # Inverse Wishart prior for the raw covariance matrix
  tau_b_raw[1:8,1:8] ~ dwish(omega_b[1:8,1:8], 9)
  sigma2_b_raw[1:8,1:8] <- inverse(tau_b_raw[,])

  ## Define elements of correlation and covariance matrices to recover
  for(k in 1:8){
    for(k_prime in 1:8){
      rho_b[k,k_prime] <- sigma2_b_raw[k,k_prime]/sqrt(sigma2_b_raw[k,k]*sigma2_b_raw[k_prime,k_prime])
      cov_b_mat[k,k_prime] <- scale_c[k]*scale_c[k_prime]*sigma2_b_raw[k,k_prime]
    }
  }

  # Variances:
  for(k in 1:8){
    var_b[k] <- cov_b_mat[k,k]
  }

  # Covariances
  cov_b[1] <- cov_b_mat[2,1]
  cov_b[2] <- cov_b_mat[3,1]
  cov_b[3] <- cov_b_mat[3,2]
  cov_b[4] <- cov_b_mat[4,1]
  cov_b[5] <- cov_b_mat[4,2]
  cov_b[6] <- cov_b_mat[4,3]
  cov_b[7] <- cov_b_mat[5,1]
  cov_b[8] <- cov_b_mat[5,2]
  cov_b[9] <- cov_b_mat[5,3]
  cov_b[10] <- cov_b_mat[5,4]
  cov_b[11] <- cov_b_mat[6,1]
  cov_b[12] <- cov_b_mat[6,2]
  cov_b[13] <- cov_b_mat[6,3]
  cov_b[14] <- cov_b_mat[6,4]
  cov_b[15] <- cov_b_mat[6,5]
  cov_b[16] <- cov_b_mat[7,1]
  cov_b[17] <- cov_b_mat[7,2]
  cov_b[18] <- cov_b_mat[7,3]
  cov_b[19] <- cov_b_mat[7,4]
  cov_b[20] <- cov_b_mat[7,5]
  cov_b[21] <- cov_b_mat[7,6]
  cov_b[22] <- cov_b_mat[8,1]
  cov_b[23] <- cov_b_mat[8,2]
  cov_b[24] <- cov_b_mat[8,3]
  cov_b[25] <- cov_b_mat[8,4]
  cov_b[26] <- cov_b_mat[8,5]
  cov_b[27] <- cov_b_mat[8,6]
  cov_b[28] <- cov_b_mat[8,7]

  # Correlations
  cor_b[1] <- rho_b[2,1]
  cor_b[2] <- rho_b[3,1]
  cor_b[3] <- rho_b[3,2]
  cor_b[4] <- rho_b[4,1]
  cor_b[5] <- rho_b[4,2]
  cor_b[6] <- rho_b[4,3]
  cor_b[7] <- rho_b[5,1]
  cor_b[8] <- rho_b[5,2]
  cor_b[9] <- rho_b[5,3]
  cor_b[10] <- rho_b[5,4]
  cor_b[11] <- rho_b[6,1]
  cor_b[12] <- rho_b[6,2]
  cor_b[13] <- rho_b[6,3]
  cor_b[14] <- rho_b[6,4]
  cor_b[15] <- rho_b[6,5]
  cor_b[16] <- rho_b[7,1]
  cor_b[17] <- rho_b[7,2]
  cor_b[18] <- rho_b[7,3]
  cor_b[19] <- rho_b[7,4]
  cor_b[20] <- rho_b[7,5]
  cor_b[21] <- rho_b[7,6]
  cor_b[22] <- rho_b[8,1]
  cor_b[23] <- rho_b[8,2]
  cor_b[24] <- rho_b[8,3]
  cor_b[25] <- rho_b[8,4]
  cor_b[26] <- rho_b[8,5]
  cor_b[27] <- rho_b[8,6]
  cor_b[28] <- rho_b[8,7]

  # Define elements of error covariance (level 1) to recover
  sigma2_error[1] <- sigma2_e[1,1]
  sigma2_error[2] <- sigma2_e[2,2]
  sigma2_error[3] <- sigma2_e[2,1]

  # Collect important parameters to assess convergence
  for_conv[1:8] <- beta_mean
  for_conv[9:16] <- var_b
  for_conv[17:44] <- cov_b
  for_conv[45:47] <- sigma2_error
  for_conv[48:75] <- cor_b
  for_conv[76] <- rho
  }"
