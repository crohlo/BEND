#' Bayesian Crossed Random Effects Model (CREM)
#'
#' @description
#' Estimates a Bayesian crossed random effects models (CREM) for longitudinal data with dynamic group membership. Four different choices for functional forms are provided: linear, quadratic, exponential, and piecewise. See Rohloff et al. (2024) for more details.
#'
#' @param data Data frame in long format, where each row describes a measurement occasion for a given individual. It is assumed that each individual has the same number of assigned timepoints (a.k.a., rows). There can be missingness in the outcome (`y_var`), but there cannot be missingness in time (`time_var`).
#' @param ind_id_var Name of column that contains ids for individuals with repeated measures in a longitudinal dataset (e.g., students).
#' @param cross_id_var Name of column that contains ids for the crossed factor (e.g., teachers).
#' @param time_var Name of column that contains the time variable. This column cannot contain any missing values.
#' @param y_var Name of column that contains the outcome variable. Missing values should be denoted by NA.
#' @param form Name of the functional form. Options include: ‘linear’ (default), ‘quadratic’, ‘exponential’, ‘piecewise’.
#' @param fixed_effects (optional) Starting values for the fixed effects parameters.
#' @param iters_adapt (optional) Number of iterations for adaptation of jags model (default = 5000).
#' @param iters_burn_in (optional) Number of iterations for burn-in (default = 50000).
#' @param iters_sampling (optional) Number of iterations for posterior sampling (default = 50000).
#' @param thin (optional) Thinning interval for posterior sampling (default = 15).
#' @param save_full_chains Logical indicating whether the MCMC chains from rjags should be saved (default = FALSE). Note, this should not be used regularly as it will result in an object with a large file size.
#' @param save_conv_chains Logical indicating whether the MCMC chains from rjags should be saved but only for the parameters monitored for convergence (default = FALSE). This would be useful for plotting traceplots for relevant model parameters to evaluate convergence behavior. Note, this should not be used regularly as it will result in an object with a large file size.
#' @param verbose Logical controlling whether progress messages/bars are generated (default = TRUE).
#'
#' @returns A list (an object of class `CREM`) with elements:
#' \item{Convergence}{Potential scale reduction factor (PSRF) for each parameter (`parameter_psrf`), Gelman multivariate scale reduction factor (`multivariate_psrf`), and mean PSRF (`mean_psrf`) to assess model convergence.}
#' \item{Model_Fit}{Deviance (`deviance`), effective number of parameters (`pD`), and Deviance information criterion (`dic`) to assess model fit.}
#' \item{Fitted_Values}{Vector giving the fitted value at each timepoint for each individual (same length as long data).}
#' \item{Functional_Form}{Functional form fitted.}
#' \item{Parameter_Estimates}{Data frame with posterior mean and 95% credible intervals for each model parameter.}
#' \item{Run_Time}{Total run time for model fitting.}
#' \item{Full_MCMC_Chains}{If save_full_chains=TRUE, raw MCMC chains from rjags.}
#' \item{Convergence_MCMC_Chains}{If save_conv_chains=TRUE, raw MCMC chains from rjags but only for the parameters monitored for convergence.}
#'
#' @details
#' For more information on the model equation and priors implemented in this function, see Rohloff et al. (2024).
#'
#' Note, this function differs from the above reference by estimating the covariances between the random effects parameters. The variance-covariance matrices of the individual and group random effects have a scaled inverse-Wishart prior (see Peralta et al., 2022).
#'
#' @author Corissa T. Rohloff
#'
#' @references
#' Peralta, Y., Kohli, N., Lock, E. F., & Davison, M. L. (2022). Bayesian modeling of associations in bivariate piecewise linear mixed-effects models. Psychological Methods, 27(1), 44–64. https://doi.org/10.1037/met0000358
#'
#' Rohloff, C. T., Kohli, N., & Lock, E. F. (2024). Identifiability and estimability of Bayesian linear and nonlinear crossed random effects models. British Journal of Mathematical and Statistical Psychology. https://doi.org/10.1111/bmsp.12334
#'
#' @examples
#' \donttest{
#' # load simulated data
#' data(SimData_PCREM)
#' # plot observed data
#' plot_BEND(data = SimData_PCREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y")
#' # fit Bayes_CREM()
#' results_pcrem <- Bayes_CREM(data = SimData_PCREM,
#'                             ind_id_var = "id",
#'                             cross_id_var = "teacherid",
#'                             time_var = "time",
#'                             y_var = "y",
#'                             form="piecewise")
#' # result summary
#' summary(results_pcrem)
#' # plot fitted results
#' plot_BEND(data = SimData_PCREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y",
#'           results = results_pcrem)
#' }
#'
#' @import stats
#'
#' @export
Bayes_CREM <- function(data,
                       ind_id_var, cross_id_var, time_var, y_var,
                       form="linear",
                       fixed_effects=NULL,
                       iters_adapt=5000, iters_burn_in=50000, iters_sampling=50000, thin=15,
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

  # Pull data for modeling ----
  y <- data[,paste0(y_var)] # vector of all outcome observations
  t <- data[,paste0(time_var)] # vector of all timepoints observations
  ind_id <- as.numeric(factor(data[,paste0(ind_id_var)])) # vector of individual ids
  cross_id <- as.numeric(factor(data[,paste0(cross_id_var)])) # vector of crossed factor ids

  # Define relevant variables ----
  n_subj <- length(unique(data[,paste0(ind_id_var)]))  # number of individuals (e.g., students)
  n_group <- length(unique(data[,paste0(cross_id_var)])) # number of groups in the crossed factor (e.g., teachers)
  n_obs <- length(y) # number of total observations

  # For piecewise model
  min_cp <- unique(sort(t))[2]
  max_cp <- unique(sort(t, decreasing = TRUE))[2]

  # Error messages for input ----

  ## data format warnings
  if(sum(length(y)==length(t))<1) stop('Columns for y_var and t_var must have the same dimensions')
  if(sum(is.na(t))>0) stop('Columns for t_var cannot have NA values (but y_var can)')
  if(min(t)!=0) warning('Prior assumes first time point measured at t=0')

  ## model specification warnings
  if(form!="linear"&&form!="quadratic"&&form!="exponential"&&form!="piecewise") stop("Specified functional form must be \"linear\", \"quadratic\", \"exponential\", or \"piecewise\"")

  # FULL MODEL ----

  ## Specify full JAGS model
  ## Linear
  if(form=="linear")      n_beta <- 2
  if(form=="linear")      full_spec <- textConnection(ln_crem)
  ## Quadratic
  if(form=="quadratic")   n_beta <- 3
  if(form=="quadratic")   full_spec <- textConnection(qd_crem)
  ## Exponential
  if(form=="exponential") n_beta <- 3
  if(form=="exponential") full_spec <- textConnection(ex_crem)
  ## Piecewise
  if(form=="piecewise")   n_beta <- 4
  if(form=="piecewise")   full_spec <- textConnection(pw_crem)

  ## Scale matrix for the Wishart distribution
  omega_b <- diag(n_beta)

  ## Variables to extract from full model
  param_recovery_full <- c('beta_ir', 'beta_mean',
                           'cov_b', 'var_b',
                           'cov_g', 'var_g',
                           'mu_y', 'sigma2_error',
                           'for_conv', 'pD', 'deviance')

  ## Create data list for JAGS model
  data_list <- list("t" = t, "y" = y,
                    "n_obs" = n_obs, "n_subj" = n_subj, 'n_group' = n_group,
                    "ind_id" = ind_id, 'cross_id' = cross_id,
                    "n_beta" = n_beta, "omega_b" = omega_b)
  if(form=="piecewise") data_list$min_cp <- min_cp
  if(form=="piecewise") data_list$max_cp <- max_cp

  ## Compile Info for Initial Values
  initial_vals <- vector('list', n_chains)
  if(!is.null(fixed_effects)){
    for(i in 1:n_chains){initial_vals[[i]]$beta_mean <- fixed_effects}
  }

  if(verbose) cat("Calibrating MCMC...\n")
  if(!is.null(fixed_effects)) full_model <- rjags::jags.model(full_spec,
                                                              data = data_list,
                                                              inits = initial_vals,
                                                              n.chains = n_chains,
                                                              n.adapt = iters_adapt,
                                                              quiet = TRUE)
  if(is.null(fixed_effects)) full_model <- rjags::jags.model(full_spec,
                                                             data = data_list,
                                                             n.chains = n_chains,
                                                             n.adapt = iters_adapt,
                                                             quiet = TRUE)

  # burn-in
  if(verbose) cat("Burn in of jags model...\n")
  update(full_model, iters_burn_in, progress.bar=progress_bar)

  # sampling
  if(verbose) cat("Collecting samples...\n")
  full_out <- rjags::jags.samples(full_model,
                                  variable.names = param_recovery_full,
                                  n.iter = iters_sampling,
                                  thin = thin,
                                  progress.bar=progress_bar)

  # Compiling Full Results -----

  if(form=="linear")      param_names <- c("beta_0_mean", "beta_1_mean",
                                           "var_b_0", "var_b_1",
                                           "cov_b_01",
                                           "var_g_0", "var_g_1",
                                           "cov_g_01",
                                           "error_var")
  if(form=="quadratic" |
     form=="exponential") param_names <- c("beta_0_mean", "beta_1_mean", "beta_2_mean",
                                           "var_b_0", "var_b_1", "var_b_2",
                                           "cov_b_01", "cov_b_02", "cov_b_12",
                                           "var_g_0", "var_g_1", "var_g_2",
                                           "cov_g_01", "cov_g_02", "cov_g_12",
                                           "error_var")

  if(form=="piecewise")   param_names <- c("beta_0_mean", "beta_1_mean", "beta_2_mean", "beta_cp_mean",
                                           "var_b_0", "var_b_1", "var_b_2", "var_b_cp",
                                           "cov_b_01", "cov_b_02", "cov_b_12", "cov_b_0cp", "cov_b_1cp", "cov_b_2cp",
                                           "var_g_0", "var_g_1", "var_g_2", "var_g_cp",
                                           "cov_g_01", "cov_g_02", "cov_g_12", "cov_g_0cp", "cov_g_1cp", "cov_g_2cp",
                                           "error_var")

  ## Convergence
  mcmc_list <- coda::as.mcmc.list(full_out$for_conv)
  for(i in 1:3){colnames(mcmc_list[[i]]) <- param_names}

  gelman_msrf <- coda::gelman.diag(mcmc_list) # individual parameter psrf, and multivariate psrf
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
  param_est <- data.frame(Mean = sum_mcmc$statistics[,1],
                          CI_Lower = sum_mcmc$quantiles[,1],
                          CI_Upper = sum_mcmc$quantiles[,5])

  ## Stop tracking run time
  run_time_total_end <- Sys.time()
  run_time_total <- run_time_total_end - run_time_total_start

  my_results <- list('Convergence' = convergence,
                     'Model_Fit' = model_fit,
                     'Fitted_Values'=y_mean,
                     'Functional_Form'=form,
                     'Parameter_Estimates'=param_est,
                     'Run_Time'=format(run_time_total))
  if(save_full_chains==TRUE){my_results$Full_MCMC_Chains=full_out}
  if(save_conv_chains==TRUE){my_results$Convergence_MCMC_Chains=mcmc_list}
  class(my_results) <- 'CREM'
  return(my_results)
}

# Models ----

## Linear ----
ln_crem <- "model{
  ##### Level-1 Model #####
  for(j in 1:n_obs){
    y[j] ~ dnorm(mu_y[j], tau_y)
    mu_y[j] <- beta_ir[j,1] + beta_ir[j,2]*t[j]
    beta_ir[j,1] <- beta_mean[1] + b0[ind_id[j],1] + g0[cross_id[j],1]
    beta_ir[j,2] <- beta_mean[2] + b0[ind_id[j],2] + g0[cross_id[j],2]
  }

  ##### Loop over Individuals #####
  for(i in 1:n_subj){
    for(p in 1:n_beta){
      # scale variances
      b0[i,p] <- b_scale_c[p]*b0_raw[i,p]
    }
    # generate raw variance-covariance matrix
    b0_raw[i,1:n_beta] ~ dmnorm(rep(0,n_beta), tau_b_raw[1:n_beta,1:n_beta])
  }

  ##### Loop over Crossed Factor #####
  for(r in 1:n_group){
    for(p in 1:n_beta){
      # scale variances
      g0[r,p] <- g_scale_c[p]*g0_raw[r,p]
    }
    # generate raw variance-covariance matrix
    g0_raw[r,1:n_beta] ~ dmnorm(rep(0,n_beta), tau_g_raw[1:n_beta,1:n_beta])
  }

  ##### Priors #####
  ## Priors for Fixed Effects ##
  for(p in 1:n_beta){
   beta_mean[p] ~ dnorm(0, 0.00001)
  }

  ## Priors for Residual Variance Components ##
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  ## Priors for Individual Random Effects ##
  tau_b_raw[1:n_beta,1:n_beta] ~ dwish(omega_b[1:n_beta,1:n_beta], n_beta+1)
  sigma2_b_raw[1:n_beta,1:n_beta] <- inverse(tau_b_raw[,])

  ## Priors for Crossed Random Effects ##
  tau_g_raw[1:n_beta,1:n_beta] ~ dwish(omega_b[1:n_beta,1:n_beta], n_beta+1)
  sigma2_g_raw[1:n_beta,1:n_beta] <- inverse(tau_g_raw[,])

  ## Priors for Scaling Constants ##
  for(p in 1:n_beta){
    # For Individual Effects #
    b_scale_c[p] ~ dunif(0, 100)
    # For Crossed Effects #
    g_scale_c[p] ~ dunif(0, 100)
  }

  ##### Variance-Covariance Matrix Elements #####
  ## Individual Effects ##
  # Scale Covariance Matrices #
  for(p in 1:n_beta){
    for(p_prime in 1:n_beta){
      cov_b_mat[p,p_prime] <- b_scale_c[p]*b_scale_c[p_prime]*sigma2_b_raw[p,p_prime]
    }
  }

  # Variances #
  for(p in 1:n_beta){
    var_b[p] <- cov_b_mat[p,p]
  }

  # Covariances #
  cov_b[1] <- cov_b_mat[2,1] # 0,1

  ## Crossed Effects ##
  # Scale Covariance Matrices #
  for(p in 1:n_beta){
    for(p_prime in 1:n_beta){
      cov_g_mat[p,p_prime] <- g_scale_c[p]*g_scale_c[p_prime]*sigma2_g_raw[p,p_prime]
    }
  }

  # Variances #
  for(p in 1:n_beta){
    var_g[p] <- cov_g_mat[p,p]
  }

  # Covariances #
  cov_g[1] <- cov_g_mat[2,1] # 0,1

  ##### Collect important parameters to assess convergence #####
  for_conv[1:2] <- beta_mean
  for_conv[3:4] <- var_b
  for_conv[5] <- cov_b
  for_conv[6:7] <- var_g
  for_conv[8] <- cov_g
  for_conv[9] <- sigma2_error
  }
"

## Quadratic ----
qd_crem <- "model{
  ##### Level-1 Model #####
  for(j in 1:n_obs){
    y[j] ~ dnorm(mu_y[j], tau_y)
    mu_y[j] <- beta_ir[j,1] + beta_ir[j,2]*t[j] + beta_ir[j,3]*(t[j]^2)
    beta_ir[j,1] <- beta_mean[1] + b0[ind_id[j],1] + g0[cross_id[j],1]
    beta_ir[j,2] <- beta_mean[2] + b0[ind_id[j],2] + g0[cross_id[j],2]
    beta_ir[j,3] <- beta_mean[3] + b0[ind_id[j],3] + g0[cross_id[j],3]
  }

  ##### Loop over Individuals #####
  for(i in 1:n_subj){
    for(p in 1:n_beta){
      # scale variances
      b0[i,p] <- b_scale_c[p]*b0_raw[i,p]
    }
    # generate raw variance-covariance matrix
    b0_raw[i,1:n_beta] ~ dmnorm(rep(0,n_beta), tau_b_raw[1:n_beta,1:n_beta])
  }

  ##### Loop over Crossed Factor #####
  for(r in 1:n_group){
    for(p in 1:n_beta){
      # scale variances
      g0[r,p] <- g_scale_c[p]*g0_raw[r,p]
    }
    # generate raw variance-covariance matrix
    g0_raw[r,1:n_beta] ~ dmnorm(rep(0,n_beta), tau_g_raw[1:n_beta,1:n_beta])
  }

  ##### Priors #####
  ## Priors for Fixed Effects ##
  for(p in 1:n_beta){
   beta_mean[p] ~ dnorm(0, 0.00001)
  }

  ## Priors for Residual Variance Components ##
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  ## Priors for Individual Random Effects ##
  tau_b_raw[1:n_beta,1:n_beta] ~ dwish(omega_b[1:n_beta,1:n_beta], n_beta+1)
  sigma2_b_raw[1:n_beta,1:n_beta] <- inverse(tau_b_raw[,])

  ## Priors for Crossed Random Effects ##
  tau_g_raw[1:n_beta,1:n_beta] ~ dwish(omega_b[1:n_beta,1:n_beta], n_beta+1)
  sigma2_g_raw[1:n_beta,1:n_beta] <- inverse(tau_g_raw[,])

  ## Priors for Scaling Constants ##
  for(p in 1:n_beta){
    # For Individual Effects #
    b_scale_c[p] ~ dunif(0, 100)
    # For Crossed Effects #
    g_scale_c[p] ~ dunif(0, 100)
  }

  ##### Variance-Covariance Matrix Elements #####
  ## Individual Effects ##
  # Scale Covariance Matrices #
  for(p in 1:n_beta){
    for(p_prime in 1:n_beta){
      cov_b_mat[p,p_prime] <- b_scale_c[p]*b_scale_c[p_prime]*sigma2_b_raw[p,p_prime]
    }
  }

  # Variances #
  for(p in 1:n_beta){
    var_b[p] <- cov_b_mat[p,p]
  }

  # Covariances #
  cov_b[1] <- cov_b_mat[2,1] # 0,1
  cov_b[2] <- cov_b_mat[3,1] # 0,2
  cov_b[3] <- cov_b_mat[3,2] # 1,2

  ## Crossed Effects ##
  # Scale Covariance Matrices #
  for(p in 1:n_beta){
    for(p_prime in 1:n_beta){
      cov_g_mat[p,p_prime] <- g_scale_c[p]*g_scale_c[p_prime]*sigma2_g_raw[p,p_prime]
    }
  }

  # Variances #
  for(p in 1:n_beta){
    var_g[p] <- cov_g_mat[p,p]
  }

  # Covariances #
  cov_g[1] <- cov_g_mat[2,1] # 0,1
  cov_g[2] <- cov_g_mat[3,1] # 0,2
  cov_g[3] <- cov_g_mat[3,2] # 1,2

  ##### Collect important parameters to assess convergence #####
  for_conv[1:3] <- beta_mean
  for_conv[4:6] <- var_b
  for_conv[7:9] <- cov_b
  for_conv[10:12] <- var_g
  for_conv[13:15] <- cov_g
  for_conv[16] <- sigma2_error
  }
"

## Exponential ----
ex_crem <- "model{
  ##### Level-1 Model #####
  for(j in 1:n_obs){
    y[j] ~ dnorm(mu_y[j], tau_y)
    mu_y[j] <- beta_ir[j,1] + beta_ir[j,2]*(1-exp(-beta_ir[j,3]*t[j]))
    beta_ir[j,1] <- beta_mean[1] + b0[ind_id[j],1] + g0[cross_id[j],1]
    beta_ir[j,2] <- beta_mean[2] + b0[ind_id[j],2] + g0[cross_id[j],2]
    beta_ir[j,3] <- beta_mean[3] + b0[ind_id[j],3] + g0[cross_id[j],3]
  }

  ##### Loop over Individuals #####
  for(i in 1:n_subj){
    for(p in 1:n_beta){
      # scale variances
      b0[i,p] <- b_scale_c[p]*b0_raw[i,p]
    }
    # generate raw variance-covariance matrix
    b0_raw[i,1:n_beta] ~ dmnorm(rep(0,n_beta), tau_b_raw[1:n_beta,1:n_beta])
  }

  ##### Loop over Crossed Factor #####
  for(r in 1:n_group){
    for(p in 1:n_beta){
      # scale variances
      g0[r,p] <- g_scale_c[p]*g0_raw[r,p]
    }
    # generate raw variance-covariance matrix
    g0_raw[r,1:n_beta] ~ dmnorm(rep(0,n_beta), tau_g_raw[1:n_beta,1:n_beta])
  }

  ##### Priors #####
  ## Priors for Fixed Effects ##
  for(p in 1:n_beta){
   beta_mean[p] ~ dnorm(0, 0.00001)
  }

  ## Priors for Residual Variance Components ##
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  ## Priors for Individual Random Effects ##
  tau_b_raw[1:n_beta,1:n_beta] ~ dwish(omega_b[1:n_beta,1:n_beta], n_beta+1)
  sigma2_b_raw[1:n_beta,1:n_beta] <- inverse(tau_b_raw[,])

  ## Priors for Crossed Random Effects ##
  tau_g_raw[1:n_beta,1:n_beta] ~ dwish(omega_b[1:n_beta,1:n_beta], n_beta+1)
  sigma2_g_raw[1:n_beta,1:n_beta] <- inverse(tau_g_raw[,])

  ## Priors for Scaling Constants ##
  for(p in 1:n_beta){
    # For Individual Effects #
    b_scale_c[p] ~ dunif(0, 100)
    # For Crossed Effects #
    g_scale_c[p] ~ dunif(0, 100)
  }

  ##### Variance-Covariance Matrix Elements #####
  ## Individual Effects ##
  # Scale Covariance Matrices #
  for(p in 1:n_beta){
    for(p_prime in 1:n_beta){
      cov_b_mat[p,p_prime] <- b_scale_c[p]*b_scale_c[p_prime]*sigma2_b_raw[p,p_prime]
    }
  }

  # Variances #
  for(p in 1:n_beta){
    var_b[p] <- cov_b_mat[p,p]
  }

  # Covariances #
  cov_b[1] <- cov_b_mat[2,1] # 0,1
  cov_b[2] <- cov_b_mat[3,1] # 0,2
  cov_b[3] <- cov_b_mat[3,2] # 1,2

  ## Crossed Effects ##
  # Scale Covariance Matrices #
  for(p in 1:n_beta){
    for(p_prime in 1:n_beta){
      cov_g_mat[p,p_prime] <- g_scale_c[p]*g_scale_c[p_prime]*sigma2_g_raw[p,p_prime]
    }
  }

  # Variances #
  for(p in 1:n_beta){
    var_g[p] <- cov_g_mat[p,p]
  }

  # Covariances #
  cov_g[1] <- cov_g_mat[2,1] # 0,1
  cov_g[2] <- cov_g_mat[3,1] # 0,2
  cov_g[3] <- cov_g_mat[3,2] # 1,2

  ##### Collect important parameters to assess convergence #####
  for_conv[1:3] <- beta_mean
  for_conv[4:6] <- var_b
  for_conv[7:9] <- cov_b
  for_conv[10:12] <- var_g
  for_conv[13:15] <- cov_g
  for_conv[16] <- sigma2_error
  }
"

## Piecewise ----
pw_crem <- "model{
  ##### Level-1 Model #####
  for(j in 1:n_obs){
    y[j] ~ dnorm(mu_y[j], tau_y)
    mu_y[j] <- beta_ir[j,1] + beta_ir[j,2]*t[j] + beta_ir[j,3]*(max(0, t[j]-beta_ir[j,4]))
    beta_ir[j,1] <- beta_mean[1] + b0[ind_id[j],1] + g0[cross_id[j],1]
    beta_ir[j,2] <- beta_mean[2] + b0[ind_id[j],2] + g0[cross_id[j],2]
    beta_ir[j,3] <- beta_mean[3] + b0[ind_id[j],3] + g0[cross_id[j],3]
    beta_ir[j,4] <- beta_mean[4] + b0[ind_id[j],4] + g0[cross_id[j],4]
  }

  ##### Loop over Individuals #####
  for(i in 1:n_subj){
    for(p in 1:n_beta){
      # scale variances
      b0[i,p] <- b_scale_c[p]*b0_raw[i,p]
    }
    # generate raw variance-covariance matrix
    b0_raw[i,1:n_beta] ~ dmnorm(rep(0,n_beta), tau_b_raw[1:n_beta,1:n_beta])
  }

  ##### Loop over Crossed Factor #####
  for(r in 1:n_group){
    for(p in 1:n_beta){
      # scale variances
      g0[r,p] <- g_scale_c[p]*g0_raw[r,p]
    }
    # generate raw variance-covariance matrix
    g0_raw[r,1:n_beta] ~ dmnorm(rep(0,n_beta), tau_g_raw[1:n_beta,1:n_beta])
  }

  ##### Priors #####
  ## Priors for Fixed Effects ##
  for(p in 1:(n_beta-1)){
   beta_mean[p] ~ dnorm(0, 0.00001)
  }
  beta_mean[n_beta] ~ dunif(min_cp,max_cp)

  ## Priors for Residual Variance Components ##
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  ## Priors for Individual Random Effects ##
  tau_b_raw[1:n_beta,1:n_beta] ~ dwish(omega_b[1:n_beta,1:n_beta], n_beta+1)
  sigma2_b_raw[1:n_beta,1:n_beta] <- inverse(tau_b_raw[,])

  ## Priors for Crossed Random Effects ##
  tau_g_raw[1:n_beta,1:n_beta] ~ dwish(omega_b[1:n_beta,1:n_beta], n_beta+1)
  sigma2_g_raw[1:n_beta,1:n_beta] <- inverse(tau_g_raw[,])

  ## Priors for Scaling Constants ##
  for(p in 1:n_beta){
    # For Individual Effects #
    b_scale_c[p] ~ dunif(0, 100)
    # For Crossed Effects #
    g_scale_c[p] ~ dunif(0, 100)
  }

  ##### Variance-Covariance Matrix Elements #####
  ## Individual Effects ##
  # Scale Covariance Matrices #
  for(p in 1:n_beta){
    for(p_prime in 1:n_beta){
      cov_b_mat[p,p_prime] <- b_scale_c[p]*b_scale_c[p_prime]*sigma2_b_raw[p,p_prime]
    }
  }

  # Variances #
  for(p in 1:n_beta){
    var_b[p] <- cov_b_mat[p,p]
  }

  # Covariances #
  cov_b[1] <- cov_b_mat[2,1] # 0,1
  cov_b[2] <- cov_b_mat[3,1] # 0,2
  cov_b[3] <- cov_b_mat[3,2] # 1,2
  cov_b[4] <- cov_b_mat[4,1] # 0,3
  cov_b[5] <- cov_b_mat[4,2] # 1,3
  cov_b[6] <- cov_b_mat[4,3] # 2,3

  ## Crossed Effects ##
  # Scale Covariance Matrices #
  for(p in 1:n_beta){
    for(p_prime in 1:n_beta){
      cov_g_mat[p,p_prime] <- g_scale_c[p]*g_scale_c[p_prime]*sigma2_g_raw[p,p_prime]
    }
  }

  # Variances #
  for(p in 1:n_beta){
    var_g[p] <- cov_g_mat[p,p]
  }

  # Covariances #
  cov_g[1] <- cov_g_mat[2,1] # 0,1
  cov_g[2] <- cov_g_mat[3,1] # 0,2
  cov_g[3] <- cov_g_mat[3,2] # 1,2
  cov_g[4] <- cov_g_mat[4,1] # 0,3
  cov_g[5] <- cov_g_mat[4,2] # 1,3
  cov_g[6] <- cov_g_mat[4,3] # 2,3

  ##### Collect important parameters to assess convergence #####
  for_conv[1:4] <- beta_mean
  for_conv[5:8] <- var_b
  for_conv[9:14] <- cov_b
  for_conv[15:18] <- var_g
  for_conv[19:24] <- cov_g
  for_conv[25] <- sigma2_error
  }
"
