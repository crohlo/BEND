#' Bayesian Piecewise Random Effects Model (PREM) + Extensions
#'
#' @description Estimates a Bayesian piecewise random effects model (PREM), with some useful extensions. There are three model options included in this function:
#' * `PREM` estimates a Bayesian piecewise random effects model with a latent number of changepoints (default). Allows the inclusion of outcome-predictive covariates (`CI-PREM`).
#' * `PREMM` estimates a piecewise random effects mixture model for a given number of latent classes and a latent number of possible changepoints in each class.
#' * `CI-PREMM` estimates a covariate influenced piecewise random effects mixture model for a given number of latent classes and a latent number of possible changepoints in each class. Allows the inclusion of outcome- and/or class-predictive covariates.
#' See Lock et al. (2018) and Lamm (2022) for more details.
#'
#' @param data Data frame in long format, where each row describes a measurement occasion for a given individual. It is assumed that each individual has the same number of assigned timepoints (a.k.a., rows). There can be missingness in the outcome (`y_var`), but there cannot be missingness in time (`time_var`).
#' @param id_var Name of column that contains ids for individuals with repeated measures in a longitudinal dataset.
#' @param time_var Name of column that contains the time variable. This column cannot contain any missing values.
#' @param y_var Name of column that contains the outcome variable. Missing values should be denoted by NA.
#' @param n_class Number of latent classes (default = 1). Note, CI-PREMM only allows for two classes.
#' @param max_cp Maximum number of changepoints in each latent class (default = 2).
#' @param class_predictive_vars Name(s) of column(s) that contain class-predictive covariates (time-invariant only). Give a vector of names if multiple covariates. Note, there cannot be any missingness in the covariates.
#' @param outcome_predictive_vars Name(s) of column(s) that contain outcome-predictive covariates (time-varying or -invariant). Give a vector of names if multiple covariates. Note, there cannot be any missingness in the covariates.
#' @param scale_prior Prior for the scale parameter for the hierarchical random effects. Options include: ‘uniform’ (scaled uniform prior; default) or ‘hc’ (scaled half-cauchy prior).
#' @param alpha Concentration parameter for Dirichlet prior for latent classes (default = 1). This can be a vector of values corresponding to the number of classes (specified by n_class). Note, this is not used for CI-PGMM.
#' @param cp_prior Prior for the number of changepoints in each class. Options include: 'binomial' (default) or 'uniform'.
#' @param binom_prob Probability for binomial prior, if specified (default = 0.5).
#' @param iters_adapt (optional) Number of iterations for adaptation of jags model (default = 1000).
#' @param iters_burn_in (optional) Number of iterations for burn-in (default = 20000).
#' @param iters_sampling (optional) Number of iterations for posterior sampling (default = 30000).
#' @param thin (optional) Thinning interval for posterior sampling (default = 15).
#' @param save_full_chains Logical indicating whether the MCMC chains from rjags should be saved (default = FALSE). Note, this should not be used regularly as it will result in an object with a large file size.
#' @param save_conv_chains Logical indicating whether the MCMC chains from rjags should be saved but only for the parameters monitored for convergence (default = FALSE). This would be useful for plotting traceplots for relevant model parameters to evaluate convergence behavior. Note, this should not be used regularly as it will result in an object with a large file size.
#' @param verbose Logical controlling whether progress messages/bars are generated (default = TRUE).
#'
#' @returns A list (an object of class `PREM`) with elements:
#' \item{Convergence}{Potential scale reduction factor (PSRF) for each parameter (`parameter_psrf`), Gelman multivariate scale reduction factor (`multivariate_psrf`), and mean PSRF (`mean_psrf`) to assess model convergence.}
#' \item{Model_Fit}{Deviance (`deviance`), effective number of parameters (`pD`), and Deviance information criterion (`dic`) to assess model fit.}
#' \item{Fitted_Values}{Vector giving the fitted value at each timepoint for each individual (same length as long data).}
#' \item{Parameter_Estimates}{Data frame with posterior mean and 95% credible intervals for each model parameter.}
#' \item{Run_Time}{Total run time for model fitting.}
#' \item{Full_MCMC_Chains}{If save_full_chains=TRUE, raw MCMC chains from rjags.}
#' \item{Convergence_MCMC_Chains}{If save_conv_chains=TRUE, raw MCMC chains from rjags but only for the parameters monitored for convergence.}
#' `Class_Information` contains a list with elements:
#' \item{class_membership}{Vector of length n with class membership assignments for each individual.}
#' \item{individ_class_probability}{nxC matrix with each individual’s probabilities of belonging to each class conditional on their class-predictive covariates (when applicable) and growth curve.}
#' \item{unconditional_class_probability}{This output will differ based on which model was fit. For a PREM or CI-PREM, this will equal 1 as there is only one class. For a PREMM or CI-PREMM with only outcome-predictive covariates, this will be a vector of length C denoting the population probability of belonging to each class. For a CI-PREMM with class-predictive covariates, this will be a vector of length n denoting the probability of each individual belonging to the non-reference class (Class 2) based on their class-predictive covariates only.}
#'
#' @details
#' For more information on the model equation and priors implemented in this function, see Lamm et al. (2022; CI-PREMM) and Lock et al. (2018; PREMM).
#'
#' @author Corissa T. Rohloff, Rik Lamm, Eric F. Lock
#'
#' @references Lamm, R. (2022). Incorporation of covariates in Bayesian piecewise growth mixture models. https://hdl.handle.net/11299/252533
#'
#' Lock, E. F., Kohli, N., & Bose, M. (2018). Detecting multiple random changepoints in Bayesian piecewise growth mixture models. Psychometrika, 83(3), 733–750. https://doi.org/10.1007/s11336-017-9594-5
#'
#' @examples
#' \donttest{
#' # load simulated data
#' data(SimData_PREM)
#' # plot observed data
#' plot_BEND(data = SimData_PREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y")
#'
#' # PREM ---------------------------------------------------------------------------------
#' # fit Bayes_PREM()
#' results_prem <- Bayes_PREM(data = SimData_PREM,
#'                            id_var = "id",
#'                            time_var = "time",
#'                            y_var = "y")
#' # result summary
#' summary(results_prem)
#' # plot fitted results
#' plot_BEND(data = SimData_PREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y",
#'           results = results_prem)
#'
#' # CI-PREM ---------------------------------------------------------------------------------
#' # fit Bayes_PREM()
#' results_ciprem <- Bayes_PREM(data = SimData_PREM,
#'                              id_var = "id",
#'                              time_var = "time",
#'                              y_var = "y",
#'                              outcome_predictive_vars = "outcome_pred_1")
#' # result summary
#' summary(results_ciprem)
#' # plot fitted results
#' plot_BEND(data = SimData_PREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y",
#'           results = results_ciprem)
#'
#' # PREMM ---------------------------------------------------------------------------------
#' # fit Bayes_PREM()
#' results_premm <- Bayes_PREM(data = SimData_PREM,
#'                             id_var = "id",
#'                             time_var = "time",
#'                             y_var = "y",
#'                             n_class = 2)
#' # result summary
#' summary(results_premm)
#' # plot fitted results
#' plot_BEND(data = SimData_PREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y",
#'           results = results_premm)
#'
#'
#' # CI-PREMM ---------------------------------------------------------------------------------
#' # fit Bayes_PREM()
#' results_cipremm <- Bayes_PREM(data = SimData_PREM,
#'                               id_var = "id",
#'                               time_var = "time",
#'                               y_var = "y",
#'                               n_class = 2,
#'                               class_predictive_vars = c("class_pred_1", "class_pred_2"),
#'                               outcome_predictive_vars = "outcome_pred_1")
#' # result summary
#' summary(results_cipremm)
#' # plot fitted results
#' plot_BEND(data = SimData_PREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y",
#'           results = results_cipremm)
#' }
#'
#' @import stats
#'
#' @export
Bayes_PREM <- function(data,
                       id_var, time_var, y_var,
                       n_class=1, max_cp=2,
                       class_predictive_vars=NULL, outcome_predictive_vars=NULL,
                       scale_prior='uniform', alpha=1,
                       cp_prior='binomial', binom_prob=0.5,
                       iters_adapt=1000, iters_burn_in=20000, iters_sampling=30000, thin=15,
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
  y <- reshape(data[,c(id_var, time_var, y_var)],
               idvar=id_var,
               timevar=time_var,
               direction='wide')
  y <- unname(as.matrix(y[,names(y)!=id_var]))

  ## time data - matrix form
  ## should be the same dimensions as y
  x <- matrix(data[,c(time_var)],
              byrow=TRUE,
              nrow=dim(y)[1],
              ncol=dim(y)[2])

  # Define relevant variables ----
  n_subj <- nrow(y)
  n_time <- ncol(y)
  max_beta <- max_cp + 2

  max_time <- max(x)
  min_time <- min(x)
  min_cp_mean <- unique(sort(x))[2]
  max_cp_mean <- unique(sort(x, decreasing=TRUE))[2]

  n_cov_class_predictive <- if(is.null(class_predictive_vars)) 0 else length(class_predictive_vars)
  n_cov_outcome_predictive <- if(is.null(outcome_predictive_vars)) 0 else length(outcome_predictive_vars)

  ## Define number of parameters (used to assess convergence later)
  n_params <- 1 + n_class*(2*max_beta) + n_class*(2*max_cp) + n_cov_outcome_predictive + n_cov_class_predictive
  if(n_cov_class_predictive>0) n_params <- n_params + 1 # only for class_predictive models (logistic intercept)

  ## for CI-PREMM only
  time_vec <- x[1,]

  # Reshape covariates (if applicable) ----

  ## class predictive covariates - matrix form
  ## assumed to be time invariant
  if(n_cov_class_predictive>0){
    class_predictive_covariates <- unique(data[,c(id_var, class_predictive_vars)])
    class_predictive_covariates <- unname(as.matrix(class_predictive_covariates[,names(class_predictive_covariates)!=id_var]))
  }

  ## outcome predictive covariates - array of matrices
  ## assumed to be time varying or invariant
  if(n_cov_outcome_predictive>0){
    outcome_predictive_covariates_list <- list()
    for(i in 1:n_cov_outcome_predictive){
      cov_matrix <- reshape(data[,c(id_var, time_var, outcome_predictive_vars[i])],
                            idvar=id_var,
                            timevar=time_var,
                            direction='wide')
      outcome_predictive_covariates_list[[i]] <- unname(as.matrix(cov_matrix[,names(cov_matrix)!=id_var]))
    }
    outcome_predictive_covariates <- simplify2array(outcome_predictive_covariates_list) # dimensions = (n_subj, n_time, n_cov_outcome_predictive)
  }

  # Error messages for input ----

  ## data format warnings
  if(sum(is.na(x)>0)) stop('Columns for t_var cannot have NA values (but y_var can)')
  if(min(x)!=0) warning('Prior assumes first time point measured at x=0')

  ## model specification warnings
  if(cp_prior!='binomial' && cp_prior!='uniform') stop('Input for cp_prior must be \'binomial\' or \'uniform\'')
  if(scale_prior!='uniform' && scale_prior!='hc') stop('Input for scale_prior must be \'uniform\' or \'hc\'')
  if(!is.null(class_predictive_vars) && !is.null(outcome_predictive_vars) && n_class>2) stop('CI-PREMM only allows for two classes')

  # Define mean & precision for random coefficients ----
  mean <- c(mean(y[,1], na.rm=TRUE), 0, rep(0,max_cp))
  prec_param_int <- 1/var(y[,1], na.rm=TRUE)
  prec_param <- (1/(max_cp*sd(y, na.rm=TRUE)/(sd(x))))^2
  prec_vec <- c(prec_param_int, rep(prec_param, max_cp+1))
  prec <- diag(max_beta)
  diag(prec) <- prec_vec

  ## Bernoulli indicator prior for uniform number of changepoints
  aux_prob <- c(max_cp:1)/c((max_cp+1):2)

  # Model specification ----
  mod_type <- if(n_class==1) ('PREM') else
    if(!is.null(class_predictive_vars) && !is.null(outcome_predictive_vars)) ('CI_PREMM_Full') else
      if(!is.null(class_predictive_vars) && is.null(outcome_predictive_vars)) ('CI_PREMM_Class_Predictive') else
        if(is.null(class_predictive_vars) && !is.null(outcome_predictive_vars)) ('CI_PREMM_Outcome_Predictive') else 'PREMM'

  # END OF SETUP ----

  if(mod_type=='PREM'){

    # PREM ----

    ## GENERATING INITIAL VALUES ----

    ## Specify fixed JAGS model
    if(cp_prior=='binomial' && n_cov_outcome_predictive==0) init_fixed_spec <- textConnection(model_prem_binomial_fixed)
    if(cp_prior=='binomial' && n_cov_outcome_predictive>0) init_fixed_spec <- textConnection(model_prem_cov_binomial_fixed)
    if(cp_prior=='uniform'  && n_cov_outcome_predictive==0) init_fixed_spec <- textConnection(model_prem_uniform_fixed)
    if(cp_prior=='uniform'  && n_cov_outcome_predictive>0) init_fixed_spec <- textConnection(model_prem_cov_uniform_fixed)

    ## Variables to extract for initialization of full model
    param_recovery_init <- c('beta', 'cp', 'n_cp', 'cp_indicator', 'sigma2_error')
    if(n_cov_outcome_predictive>0) param_recovery_init <- c(param_recovery_init, 'outcome_predictive_covariate_alpha')

    ## Create data list for JAGS model
    data_list <- list('x'=x, 'y'=y,
                      'n_subj'=n_subj, 'n_time'=n_time,
                      'max_cp'=max_cp, 'max_beta'=max_beta,
                      'mean'=mean, 'prec'=prec,
                      'min_cp_mean'=min_cp_mean, 'max_cp_mean'=max_cp_mean,
                      'n_class'=n_class)
    if(cp_prior=='binomial') data_list$binom_prob <- binom_prob
    if(cp_prior=='uniform') data_list$aux_prob <- aux_prob
    if(n_cov_outcome_predictive>0) data_list$n_cov_outcome_predictive <- n_cov_outcome_predictive
    if(n_cov_outcome_predictive>0) data_list$outcome_predictive_covariates <- outcome_predictive_covariates

    if(verbose) cat('Computing initial values...\n')
    init_fixed_model <- rjags::jags.model(init_fixed_spec,
                                          data=data_list,
                                          n.chains=n_chains,
                                          n.adapt=500,
                                          quiet=TRUE)
    ## burn-in
    update(init_fixed_model, 2000, progress.bar=progress_bar)

    ## sampling
    init_fixed_out <- rjags::jags.samples(init_fixed_model,
                                          variable.names=param_recovery_init,
                                          n.iter=2000,
                                          thin=4,
                                          progress.bar=progress_bar)

    ### Compiling Initialization Results ----

    ## Permute changepoint labels, if necessary, so they are ordered correctly
    init_fixed_out <- Permute_CP(output=init_fixed_out,
                                 n_class=n_class,
                                 max_cp=max_cp,
                                 max_beta=max_beta)

    ## Extract initial values for full model parameters
    initial_vals <- vector('list', n_chains)
    for(i in 1:n_chains){initial_vals[[i]]$beta_mean <- matrix(nrow=n_class, ncol=max_beta)}
    for(i in 1:n_chains){initial_vals[[i]]$beta_sd <- matrix(nrow=n_class, ncol=max_beta)}
    for(i in 1:n_chains){initial_vals[[i]]$cp_mean <- matrix(nrow=n_class, ncol=max_cp)}
    for(i in 1:n_chains){initial_vals[[i]]$cp_sd <- matrix(nrow=n_class, ncol=max_cp)}

    for(i in 1:n_chains){
      for(j in 1:n_class){
        initial_vals[[i]]$beta_mean[j,] <- rowMeans(init_fixed_out$beta[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
        if(length(initial_vals[[i]]$cp_mean)<2){
          initial_vals[[i]]$cp_mean[j,] <- mean(init_fixed_out$cp[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
        }
        if(length(initial_vals[[i]]$cp_mean)>1){
          initial_vals[[i]]$cp_mean[j,] <- rowMeans(init_fixed_out$cp[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
        }
        for(k in 1:max_cp){
          initial_vals[[i]]$cp_sd[j,k] <- (max_time-min_time)/(4*max_cp)
        }
        initial_vals[[i]]$beta_sd[j,1] <- sqrt(1/prec_param_int)
        for(k in 2:max_beta){
          initial_vals[[i]]$beta_sd[j,k] <- sqrt(2/prec_param)
        }
      }
    }

    ## only when cp_prior=='binomial'
    if(cp_prior=='binomial'){
      for(i in 1:n_chains){initial_vals[[i]]$cp_indicator <- matrix(nrow=n_class, ncol=max_cp)}
      for(i in 1:n_chains){
        for(j in 1:n_class){
          for(k in 1:max_cp){
            initial_vals[[i]]$cp_indicator[j,k] <- getmode(init_fixed_out$cp_indicator[j,k,,i])
          }
        }
      }
    }

    ## only when cp_prior=='uniform'
    if(cp_prior=='uniform'){
      for(i in 1:n_chains){initial_vals[[i]]$Temp <- matrix(nrow=n_class, ncol=max_cp)}
      for(i in 1:n_chains){
        for(j in 1:n_class){
          for(k in 1:max_cp){
            initial_vals[[i]]$Temp[j,k] <- getmode(init_fixed_out$cp_indicator[j,k,,i])
          }
        }
      }
    }

    ## only when n_cov_outcome_predictive>0
    if(n_cov_outcome_predictive>0){
      for(i in 1:n_chains){
        if(n_cov_outcome_predictive==1){
          initial_vals[[i]]$outcome_predictive_covariate_alpha <- mean(init_fixed_out$outcome_predictive_covariate_alpha[,,i])
        }
        if(n_cov_outcome_predictive>1){
          initial_vals[[i]]$outcome_predictive_covariate_alpha <- rowMeans(init_fixed_out$outcome_predictive_covariate_alpha[,,i])
        }
      }
    }

    ## FULL MODEL -----

    ## Specify full JAGS model
    # without outcome predictive covariates
    if(n_cov_outcome_predictive==0 && cp_prior=='binomial' && scale_prior=='uniform') full_spec <- textConnection(model_prem_binomial_scaleunif)
    if(n_cov_outcome_predictive==0 && cp_prior=='binomial' && scale_prior=='hc') full_spec <- textConnection(model_prem_binomial_scalehc)
    if(n_cov_outcome_predictive==0 && cp_prior=='uniform' && scale_prior=='uniform') full_spec <- textConnection(model_prem_uniform_scaleunif)
    if(n_cov_outcome_predictive==0 && cp_prior=='uniform' && scale_prior=='hc') full_spec <- textConnection(model_prem_uniform_scalehc)
    # with outcome predictive covariates
    if(n_cov_outcome_predictive>0 && cp_prior=='binomial' && scale_prior=='uniform') full_spec <- textConnection(model_prem_cov_binomial_scaleunif)
    if(n_cov_outcome_predictive>0 && cp_prior=='binomial' && scale_prior=='hc') full_spec <- textConnection(model_prem_cov_binomial_scalehc)
    if(n_cov_outcome_predictive>0 && cp_prior=='uniform' && scale_prior=='uniform') full_spec <- textConnection(model_prem_cov_uniform_scaleunif)
    if(n_cov_outcome_predictive>0 && cp_prior=='uniform' && scale_prior=='hc') full_spec <- textConnection(model_prem_cov_uniform_scalehc)

    ## Variables to extract from full model
    param_recovery_full <- c('beta', 'beta_mean', 'beta_sd',
                             'cp', 'cp_mean', 'cp_sd',
                             'n_cp', 'cp_indicator',
                             'mu_y', 'sigma2_error',
                             'for_conv','pD','deviance')
    if(n_cov_outcome_predictive>0) param_recovery_full <- c(param_recovery_full, 'outcome_predictive_covariate_alpha')

    ## Create data list for JAGS model
    data_list <- list('x'=x, 'y'=y,
                      'n_subj'=n_subj, 'n_time'=n_time,
                      'min_time'=min_time, 'max_time'=max_time,
                      'max_cp'=max_cp, 'max_beta'=max_beta,
                      'mean'=mean, 'prec'=prec,
                      'min_cp_mean'=min_cp_mean, 'max_cp_mean'=max_cp_mean,
                      'n_class'=n_class,
                      'n_params'=n_params)
    if(cp_prior=='binomial') data_list$binom_prob <- binom_prob
    if(cp_prior=='uniform') data_list$aux_prob <- aux_prob
    if(n_cov_outcome_predictive>0) data_list$n_cov_outcome_predictive <- n_cov_outcome_predictive
    if(n_cov_outcome_predictive>0) data_list$outcome_predictive_covariates <- outcome_predictive_covariates

    if(verbose) cat('Calibrating MCMC...\n')
    full_model <- rjags::jags.model(full_spec,
                                    data=data_list,
                                    inits=initial_vals,
                                    n.chains=n_chains,
                                    n.adapt=iters_adapt,
                                    quiet=TRUE)

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

    ## Permute changepoint labels, if necessary, so they are ordered correctly
    full_out <- Permute_CP(output=full_out,
                           n_class=n_class,
                           max_cp=max_cp,
                           max_beta=max_beta)
    #####

  } else
    if(mod_type=='PREMM'){

      # PREMM ----

      ## GENERATING INITIAL VALUES ----

      ## Specify fixed JAGS model
      if(cp_prior=='binomial') init_fixed_spec <- textConnection(model_premm_binomial_fixed)
      if(cp_prior=='uniform')  init_fixed_spec <- textConnection(model_premm_uniform_fixed)

      ## Variables to extract for initialization of full model
      param_recovery_init <- c('beta', 'cp', 'n_cp', 'cp_indicator', 'sigma2_error', 'class')

      ## Create data list for JAGS model
      data_list <- list('x'=x, 'y'=y,
                        'n_subj'=n_subj, 'n_time'=n_time,
                        'max_cp'=max_cp, 'max_beta'=max_beta,
                        'mean'=mean, 'prec'=prec,
                        'min_cp_mean'=min_cp_mean, 'max_cp_mean'=max_cp_mean,
                        'n_class'=n_class, 'alpha'=rep(alpha, n_class))
      if(cp_prior=='binomial') data_list$binom_prob <- binom_prob
      if(cp_prior=='uniform') data_list$aux_prob <- aux_prob

      if(verbose) cat('Computing initial values...\n')
      init_fixed_model <- rjags::jags.model(init_fixed_spec,
                                            data=data_list,
                                            n.chains=n_chains,
                                            n.adapt=500,
                                            quiet=TRUE)
      ## burn-in
      update(init_fixed_model, 2000, progress.bar=progress_bar)

      ## sampling
      init_fixed_out <- rjags::jags.samples(init_fixed_model,
                                            variable.names=param_recovery_init,
                                            n.iter=2000,
                                            thin=4,
                                            progress.bar=progress_bar)

      ### Compiling Initialization Results ----

      ## Permute changepoint labels, if necessary, so they are ordered correctly
      init_fixed_out <- Permute_CP(output=init_fixed_out,
                                   n_class=n_class,
                                   max_cp=max_cp,
                                   max_beta=max_beta)

      ## Realign Classes - checking if the classes are properly ordered
      init_fixed_out <- Realign_ECR(output=init_fixed_out,
                                    n_class=n_class,
                                    model=mod_type)

      ## Extract initial values for full model parameters
      initial_vals <- vector('list', n_chains)
      for(i in 1:n_chains){initial_vals[[i]]$beta_mean <- matrix(nrow=n_class, ncol=max_beta)}
      for(i in 1:n_chains){initial_vals[[i]]$beta_sd <- matrix(nrow=n_class, ncol=max_beta)}
      for(i in 1:n_chains){initial_vals[[i]]$cp_mean <- matrix(nrow=n_class, ncol=max_cp)}
      for(i in 1:n_chains){initial_vals[[i]]$cp_sd <- matrix(nrow=n_class, ncol=max_cp)}
      for(i in 1:n_chains){initial_vals[[i]]$class <- NA}

      for(i in 1:n_chains){
        initial_vals[[i]]$class <- apply(init_fixed_out$class[,,i],1,getmode)
        for(j in 1:n_class){
          initial_vals[[i]]$beta_mean[j,] <- rowMeans(init_fixed_out$beta[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
          if(length(initial_vals[[i]]$cp_mean)<2){
            initial_vals[[i]]$cp_mean[j,] <- mean(init_fixed_out$cp[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
          }
          if(length(initial_vals[[i]]$cp_mean)>1){
            initial_vals[[i]]$cp_mean[j,] <- rowMeans(init_fixed_out$cp[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
          }
          for(k in 1:max_cp){
            initial_vals[[i]]$cp_sd[j,k] <- (max_time-min_time)/(4*max_cp)
          }
          initial_vals[[i]]$beta_sd[j,1] <- sqrt(1/prec_param_int)
          for(k in 2:max_beta){
            initial_vals[[i]]$beta_sd[j,k] <- sqrt(2/prec_param)
          }
        }
      }

      ## only when cp_prior=='binomial'
      if(cp_prior=='binomial'){
        for(i in 1:n_chains){initial_vals[[i]]$cp_indicator <- matrix(nrow=n_class, ncol=max_cp)}
        for(i in 1:n_chains){
          for(j in 1:n_class){
            for(k in 1:max_cp){
              initial_vals[[i]]$cp_indicator[j,k] <- getmode(init_fixed_out$cp_indicator[j,k,,i])
            }
          }
        }
      }

      ## only when cp_prior=='uniform'
      if(cp_prior=='uniform'){
        for(i in 1:n_chains){initial_vals[[i]]$Temp <- matrix(nrow=n_class, ncol=max_cp)}
        for(i in 1:n_chains){
          for(j in 1:n_class){
            for(k in 1:max_cp){
              initial_vals[[i]]$Temp[j,k] <- getmode(init_fixed_out$cp_indicator[j,k,,i])
            }
          }
        }
      }

      ## FULL MODEL -----

      ## Specify full JAGS model
      if(cp_prior=='binomial' && scale_prior=='uniform') full_spec <- textConnection(model_premm_binomial_scaleunif)
      if(cp_prior=='binomial' && scale_prior=='hc') full_spec <- textConnection(model_premm_binomial_scalehc)
      if(cp_prior=='uniform' && scale_prior=='uniform') full_spec <- textConnection(model_premm_uniform_scaleunif)
      if(cp_prior=='uniform' && scale_prior=='hc') full_spec <- textConnection(model_premm_uniform_scalehc)

      ## Variables to extract from full model
      param_recovery_full <- c('beta', 'beta_mean', 'beta_sd',
                               'cp', 'cp_mean', 'cp_sd',
                               'n_cp', 'cp_indicator',
                               'class_prob', 'class',
                               'mu_y', 'sigma2_error',
                               'for_conv','pD','deviance')

      ## Create data list for JAGS model
      data_list <- list('x'=x, 'y'=y,
                        'n_subj'=n_subj, 'n_time'=n_time,
                        'min_time'=min_time, 'max_time'=max_time,
                        'max_cp'=max_cp, 'max_beta'=max_beta,
                        'mean'=mean, 'prec'=prec,
                        'min_cp_mean'=min_cp_mean, 'max_cp_mean'=max_cp_mean,
                        'n_class'=n_class, 'alpha'=rep(alpha, n_class),
                        'n_params'=n_params)
      if(cp_prior=='binomial') data_list$binom_prob <- binom_prob
      if(cp_prior=='uniform') data_list$aux_prob <- aux_prob

      if(verbose) cat('Calibrating MCMC...\n')
      full_model <- rjags::jags.model(full_spec,
                                      data=data_list,
                                      inits=initial_vals,
                                      n.chains=n_chains,
                                      n.adapt=iters_adapt,
                                      quiet=TRUE)

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

      ## Permute changepoint labels, if necessary, so they are ordered correctly
      full_out <- Permute_CP(output=full_out,
                             n_class=n_class,
                             max_cp=max_cp,
                             max_beta=max_beta)

      ## Realign Classes - checking if the classes are properly ordered
      full_out <- Realign_ECR(output=full_out,
                              n_class=n_class,
                              model=mod_type)
      #####

    } else
    {

      # CI-PREMM ----

      ## GENERATE INITIAL VALUES ----

      ## Specify fixed JAGS model
      if(mod_type=='CI_PREMM_Full' && cp_prior=='binomial') init_fixed_spec <- textConnection(model_cipremm_full_binomial_fixed)
      if(mod_type=='CI_PREMM_Full' && cp_prior=='uniform') init_fixed_spec <- textConnection(model_cipremm_full_uniform_fixed)
      if(mod_type=='CI_PREMM_Class_Predictive' && cp_prior=='binomial') init_fixed_spec <- textConnection(model_cipremm_cpo_binomial_fixed)
      if(mod_type=='CI_PREMM_Class_Predictive' && cp_prior=='uniform') init_fixed_spec <- textConnection(model_cipremm_cpo_uniform_fixed)
      if(mod_type=='CI_PREMM_Outcome_Predictive' && cp_prior=='binomial') init_fixed_spec <- textConnection(model_cipremm_opo_binomial_fixed)
      if(mod_type=='CI_PREMM_Outcome_Predictive' && cp_prior=='uniform') init_fixed_spec <- textConnection(model_cipremm_opo_uniform_fixed)

      ## Variables to extract for initialization of full model
      param_recovery_init <- c('beta', 'cp', 'n_cp', 'cp_indicator', 'sigma2_error', 'class')
      if(mod_type=='CI_PREMM_Full') param_recovery_init <- c(param_recovery_init, c('logistic_intercept', 'class_predictive_covariate_lambda', 'outcome_predictive_covariate_alpha'))
      if(mod_type=='CI_PREMM_Class_Predictive') param_recovery_init <- c(param_recovery_init, c('logistic_intercept', 'class_predictive_covariate_lambda'))
      if(mod_type=='CI_PREMM_Outcome_Predictive') param_recovery_init <- c(param_recovery_init, 'outcome_predictive_covariate_alpha')

      ## Create data list for JAGS model
      data_list <- list('x'=x, 'y'=y,
                        'n_subj'=n_subj, 'n_time'=n_time,
                        'max_cp'=max_cp, 'max_beta'=max_beta,
                        'mean'=mean, 'prec'=prec,
                        'min_cp_mean'=min_cp_mean, 'max_cp_mean'=max_cp_mean,
                        'n_class'=n_class)
      if(cp_prior=='binomial') data_list$binom_prob <- binom_prob
      if(cp_prior=='uniform') data_list$aux_prob <- aux_prob
      if(mod_type=='CI_PREMM_Outcome_Predictive') data_list$alpha <- rep(alpha, n_class)
      if(mod_type=='CI_PREMM_Full' | mod_type=='CI_PREMM_Class_Predictive') data_list$class_predictive_covariates <- class_predictive_covariates
      if(mod_type=='CI_PREMM_Full' | mod_type=='CI_PREMM_Class_Predictive') data_list$n_cov_class_predictive <- n_cov_class_predictive
      if(mod_type=='CI_PREMM_Full' | mod_type=='CI_PREMM_Outcome_Predictive') data_list$outcome_predictive_covariates <- outcome_predictive_covariates
      if(mod_type=='CI_PREMM_Full' | mod_type=='CI_PREMM_Outcome_Predictive') data_list$n_cov_outcome_predictive <- n_cov_outcome_predictive

      if(verbose) cat('Computing initial values...\n')
      init_fixed_model <- rjags::jags.model(init_fixed_spec,
                                            data=data_list,
                                            n.chains=n_chains,
                                            n.adapt=500,
                                            quiet=TRUE)
      ## burn-in
      update(init_fixed_model, 2000, progress.bar=progress_bar)

      ## sampling
      init_fixed_out <- rjags::jags.samples(init_fixed_model,
                                            variable.names=param_recovery_init,
                                            n.iter=2000,
                                            thin=4,
                                            progress.bar=progress_bar)

      ### Compiling Initialization Results ----

      ## Permute changepoint labels, if necessary, so they are ordered correctly
      init_fixed_out <- Permute_CP(output=init_fixed_out,
                                   n_class=n_class,
                                   max_cp=max_cp,
                                   max_beta=max_beta)

      ## Realign Classes - checking if the classes are properly ordered
      init_fixed_out <- Realign_ECR(output=init_fixed_out,
                                    n_class=n_class,
                                    model=mod_type)

      ## Extract initial values for full model parameters
      initial_vals <- vector('list', n_chains)
      for(i in 1:n_chains){initial_vals[[i]]$beta_mean <- matrix(nrow=n_class, ncol=max_beta)}
      for(i in 1:n_chains){initial_vals[[i]]$beta_sd <- matrix(nrow=n_class, ncol=max_beta)}
      for(i in 1:n_chains){initial_vals[[i]]$cp_mean <- matrix(nrow=n_class, ncol=max_cp)}
      for(i in 1:n_chains){initial_vals[[i]]$cp_sd <- matrix(nrow=n_class, ncol=max_cp)}

      for(i in 1:n_chains){
        for(j in 1:n_class){
          initial_vals[[i]]$beta_mean[j,] <- rowMeans(init_fixed_out$beta[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
          if(length(initial_vals[[i]]$cp_mean)<2){
            initial_vals[[i]]$cp_mean[j,] <- mean(init_fixed_out$cp[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
          }
          if(length(initial_vals[[i]]$cp_mean)>1){
            initial_vals[[i]]$cp_mean[j,] <- rowMeans(init_fixed_out$cp[j,,init_fixed_out$n_cp[j,,i]==getmode(init_fixed_out$n_cp[j,,i]),i])
          }
          for(k in 1:max_cp){
            initial_vals[[i]]$cp_sd[j,k] <- (max_time-min_time)/(4*max_cp)
          }
          initial_vals[[i]]$beta_sd[j,1] <- sqrt(1/prec_param_int)
          for(k in 2:max_beta){
            initial_vals[[i]]$beta_sd[j,k] <- sqrt(2/prec_param)
          }
        }
      }

      ## only for CI_PREMM_Outcome_Predictive
      if(mod_type=='CI_PREMM_Outcome_Predictive'){
        for(i in 1:n_chains){
          initial_vals[[i]]$class <- apply(init_fixed_out$class[,,i],1,getmode)
        }
      }

      ## only when cp_prior=='binomial'
      if(cp_prior=='binomial'){
        for(i in 1:n_chains){initial_vals[[i]]$cp_indicator <- matrix(nrow=n_class, ncol=max_cp)}
        for(i in 1:n_chains){
          for(j in 1:n_class){
            for(k in 1:max_cp){
              initial_vals[[i]]$cp_indicator[j,k] <- getmode(init_fixed_out$cp_indicator[j,k,,i])
            }
          }
        }
      }

      ## only when cp_prior=='uniform'
      if(cp_prior=='uniform'){
        for(i in 1:n_chains){initial_vals[[i]]$Temp <- matrix(nrow=n_class, ncol=max_cp)}
        for(i in 1:n_chains){
          for(j in 1:n_class){
            for(k in 1:max_cp){
              initial_vals[[i]]$Temp[j,k] <- getmode(init_fixed_out$cp_indicator[j,k,,i])
            }
          }
        }
      }

      ## only when n_cov_outcome_predictive>0
      if(n_cov_outcome_predictive>0){
        for(i in 1:n_chains){
          if(n_cov_outcome_predictive==1){
            initial_vals[[i]]$outcome_predictive_covariate_alpha <- mean(init_fixed_out$outcome_predictive_covariate_alpha[,,i])
          }
          if(n_cov_outcome_predictive>1){
            initial_vals[[i]]$outcome_predictive_covariate_alpha <- rowMeans(init_fixed_out$outcome_predictive_covariate_alpha[,,i])
          }
        }
      }

      ## only when n_cov_class_predictive>0
      if(n_cov_class_predictive>0){
        for(i in 1:n_chains){
          if(n_cov_class_predictive==1){
            initial_vals[[i]]$class_predictive_covariate_lambda <- mean(init_fixed_out$class_predictive_covariate_lambda[,,i])
          }
          if(n_cov_class_predictive>1){
            initial_vals[[i]]$class_predictive_covariate_lambda <- rowMeans(init_fixed_out$class_predictive_covariate_lambda[,,i])
          }
          initial_vals[[i]]$logistic_intercept <- mean(init_fixed_out$logistic_intercept[,,i])
        }
      }

      ## FULL MODEL -----

      ## Specify full JAGS model
      if(mod_type=='CI_PREMM_Full') model_list <- cipremm_mods_full
      if(mod_type=='CI_PREMM_Class_Predictive') model_list <- cipremm_mods_cpo
      if(mod_type=='CI_PREMM_Outcome_Predictive') model_list <- cipremm_mods_opo

      if(cp_prior=='binomial' && scale_prior=='uniform') full_spec <- textConnection(model_list[[1]])
      if(cp_prior=='binomial' && scale_prior=='hc') full_spec <- textConnection(model_list[[2]])
      if(cp_prior=='uniform' && scale_prior=='uniform') full_spec <- textConnection(model_list[[3]])
      if(cp_prior=='uniform' && scale_prior=='hc') full_spec <- textConnection(model_list[[4]])

      ## Variables to extract from full model
      param_recovery_full <- c('beta', 'beta_mean', 'beta_sd',
                               'cp', 'cp_mean', 'cp_sd',
                               'n_cp', 'cp_indicator',
                               'class',
                               'mu_y', 'sigma2_error',
                               'for_conv','pD','deviance')
      if(mod_type=='CI_PREMM_Full') param_recovery_full <- c(param_recovery_full, c('logistic_intercept', 'class_predictive_covariate_lambda', 'outcome_predictive_covariate_alpha', 'logistic_class_prob'))
      if(mod_type=='CI_PREMM_Class_Predictive') param_recovery_full <- c(param_recovery_full, c('logistic_intercept', 'class_predictive_covariate_lambda', 'logistic_class_prob'))
      if(mod_type=='CI_PREMM_Outcome_Predictive') param_recovery_full <- c(param_recovery_full, 'outcome_predictive_covariate_alpha', 'class_prob')

      ## Create data list for JAGS model
      data_list <- list('x'=x, 'y'=y,
                        'n_subj'=n_subj, 'n_time'=n_time,
                        'min_time'=min_time, 'max_time'=max_time,
                        'max_cp'=max_cp, 'max_beta'=max_beta,
                        'mean'=mean, 'prec'=prec,
                        'min_cp_mean'=min_cp_mean, 'max_cp_mean'=max_cp_mean,
                        'n_class'=n_class,
                        'n_params'=n_params)
      if(cp_prior=='binomial') data_list$binom_prob <- binom_prob
      if(cp_prior=='uniform') data_list$aux_prob <- aux_prob
      if(mod_type=='CI_PREMM_Outcome_Predictive') data_list$alpha <- rep(alpha, n_class)
      if(mod_type=='CI_PREMM_Full' | mod_type=='CI_PREMM_Class_Predictive') data_list$class_predictive_covariates <- class_predictive_covariates
      if(mod_type=='CI_PREMM_Full' | mod_type=='CI_PREMM_Class_Predictive') data_list$n_cov_class_predictive <- n_cov_class_predictive
      if(mod_type=='CI_PREMM_Full' | mod_type=='CI_PREMM_Outcome_Predictive') data_list$outcome_predictive_covariates <- outcome_predictive_covariates
      if(mod_type=='CI_PREMM_Full' | mod_type=='CI_PREMM_Outcome_Predictive') data_list$n_cov_outcome_predictive <- n_cov_outcome_predictive

      if(verbose) cat('Calibrating MCMC...\n')
      full_model <- rjags::jags.model(full_spec,
                                      data=data_list,
                                      inits=initial_vals,
                                      n.chains=n_chains,
                                      n.adapt=iters_adapt,
                                      quiet=TRUE)

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

      ## Permute changepoint labels, if necessary, so they are ordered correctly
      full_out <- Permute_CP(output=full_out,
                             n_class=n_class,
                             max_cp=max_cp,
                             max_beta=max_beta)

      ## Realign Classes - checking if the classes are properly ordered
      full_out <- Realign_ECR(output=full_out,
                              n_class=n_class,
                              model=mod_type)
    }

  # Compiling Full Results -----

  ## Convergence
  # growth parameters (class dependent)
  for(c in 1:n_class){
    full_out$for_conv[1:max_beta+(4+4*max_cp)*(c-1),,] <- full_out$beta_mean[c,,,]
    full_out$for_conv[(1+max_beta):(2*max_beta)+(4+4*max_cp)*(c-1),,] <- full_out$beta_sd[c,,,]
    full_out$for_conv[(1+2*max_beta):(2*max_beta+max_cp)+(4+4*max_cp)*(c-1),,] <- full_out$cp_mean[c,,,]
    full_out$for_conv[(1+2*max_beta+max_cp):(2*max_beta+2*max_cp)+(4+4*max_cp)*(c-1),,] <- full_out$cp_sd[c,,,]
  }
  # covariates (if applicable)
  if(n_cov_outcome_predictive>0){
    full_out$for_conv[(1+2*max_beta+2*max_cp):(2*max_beta+2*max_cp+n_cov_outcome_predictive)+(4+4*max_cp)*(n_class-1),,] <- full_out$outcome_predictive_covariate_alpha[,,]
  }
  if(n_cov_class_predictive>0){
    full_out$for_conv[(1+2*max_beta+2*max_cp+n_cov_outcome_predictive):(2*max_beta+2*max_cp+n_cov_outcome_predictive+n_cov_class_predictive)+(4+4*max_cp)*(n_class-1),,] <- full_out$class_predictive_covariate_lambda[,,]
  }
  if(n_cov_class_predictive>0){
    full_out$for_conv[(1+2*max_beta+2*max_cp+n_cov_outcome_predictive+n_cov_class_predictive)+(4+4*max_cp)*(n_class-1),,] <- full_out$logistic_intercept[,,]
  }
  # error variance
  full_out$for_conv[n_params,,] <- full_out$sigma2_error

  mcmc_list <- coda::as.mcmc.list(full_out$for_conv)
  gelman_msrf <- coda::gelman.diag(mcmc_list)

  # Initialize the parameter vector
  param_names <- c('NA')
  for(c in 1:n_class){
    param_names <- c(param_names,
                     paste0('Class ', c, ': beta_', 0:(max_beta-1),'_mean'),
                     paste0('Class ', c, ': beta_', 0:(max_beta-1),'_sd'),
                     paste0('Class ', c, ': cp_', 1:max_cp,'_mean'),
                     paste0('Class ', c, ': cp_', 1:max_cp,'_sd'))
  }
  param_names <- c(param_names[-1],
                   if(n_cov_outcome_predictive>0) paste0('outcome_predictive_covariate_alpha_', 1:n_cov_outcome_predictive),
                   if(n_cov_class_predictive>0) paste0('class_predictive_covariate_lambda_', 1:n_cov_class_predictive),
                   if(n_cov_class_predictive>0) 'logistic_intercept',
                   'error_var')
  row.names(gelman_msrf$psrf) <- param_names

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

  ## Class Information
  # each of these is fixed to 1 for PREM
  class_membership <- rep(1, n_subj)
  individ_class_probability <- matrix(1, nrow=n_subj, ncol=n_class) # conditional on logistic model (when applicable) and growth curve
  unconditional_class_probability <- 1
  if(n_class>1){
    for(i in 1:n_subj){class_membership[i] <- as.numeric(names(which.max(table(full_out$class[i,,]))))}
    for(i in 1:n_subj){individ_class_probability[i,1:n_class] <- prop.table(table(factor(full_out$class[i,,], levels = 1:n_class)))}
    if(n_cov_class_predictive==0){
      unconditional_class_probability <- apply(full_out$class_prob, c(1), FUN='mean') # not based on growth curve part of the model
    }
    if(n_cov_class_predictive>0){
      unconditional_class_probability <- apply(full_out$logistic_class_prob, c(1), FUN='mean')
    }
  }
  class_info <- list('class_membership'=class_membership,
                     'individ_class_probability'=individ_class_probability,
                     'unconditional_class_probability'=unconditional_class_probability)

  ## Parameter Estimates
  param_est <- list()
  for(c in 1:n_class){

    # Probability for number of changepoints
    K_prob <- matrix(nrow=(max_cp+1), ncol=1)
    for(k in 0:max_cp){
      K_prob[k+1] <- sum(full_out$n_cp[c,,]==k)/length(full_out$n_cp[c,,])
    }
    colnames(K_prob) <- 'Probability'
    rownames(K_prob) <- c(paste0('K=',0:max_cp))

    # Parameter estimates for each result for the number of changepoints
    K=list()
    for(k in 0:max_cp){
      K[[k+1]] <- list()
      if(K_prob[k+1]>0.01){
        beta_array <- array(full_out$beta, dim=c(dim(full_out$beta)[1],dim(full_out$beta)[2],dim(full_out$beta)[3]*dim(full_out$beta)[4]))
        beta_array <- beta_array[,,full_out$n_cp[c,,]==k]
        beta <- apply(beta_array, c(1,2), 'mean')[,1:(2+k)]
        beta_mean_array <- array(full_out$beta_mean[c,,,], dim=c(dim(full_out$beta_mean)[2],dim(full_out$beta_mean)[3]*dim(full_out$beta_mean)[4]))
        beta_mean_array <- beta_mean_array[1:(2+k),full_out$n_cp[c,,]==k]
        b_mean <- apply(beta_mean_array, c(1), 'mean')
        b_mean_CI <- apply(beta_mean_array,c(1), 'quantile', probs=c(0.025,0.975))
        beta_sd_array <- array(full_out$beta_sd[c,,,], dim=c(dim(full_out$beta_sd)[2],dim(full_out$beta_sd)[3]*dim(full_out$beta_sd)[4]))
        beta_sd_array <- beta_sd_array[1:(2+k),full_out$n_cp[c,,]==k]
        beta_sd <- apply(beta_sd_array, c(1), 'mean')
        beta_sd_CI <- apply(beta_sd_array, c(1), 'quantile', probs=c(0.025,0.975))
        cp <- cp_mean <- cp_mean_CI <- cp_sd <- cp_sd_CI <- NULL
        if(k>0){
          cp_array <- array(full_out$cp, dim=c(dim(full_out$cp)[1],dim(full_out$cp)[2],dim(full_out$cp)[3]*dim(full_out$cp)[4]))
          cp_array <- cp_array[,,full_out$n_cp[c,,]==k]
          cp <- apply(cp_array,c(1,2),'mean')[,1:k]
          cp_mean_array <- array(full_out$cp_mean[c,,,], dim=c(dim(full_out$cp_mean)[2],dim(full_out$cp_mean)[3]*dim(full_out$cp_mean)[4]))
          cp_mean_array <- cp_mean_array[1:k,full_out$n_cp[c,,]==k,drop=FALSE]
          cp_mean <- apply(cp_mean_array,c(1),'mean')
          cp_mean_CI <- apply(cp_mean_array,c(1),'quantile',probs=c(0.025,0.975))
          cp_sd_array <- array(full_out$cp_sd[c,,,], dim=c(dim(full_out$cp_sd)[2],dim(full_out$cp_sd)[3]*dim(full_out$cp_sd)[4]))
          cp_sd_array <- cp_sd_array[1:k,full_out$n_cp[c,,]==k,drop=FALSE]
          cp_sd <- apply(cp_sd_array,c(1),'mean')
          cp_sd_CI <- apply(cp_sd_array,c(1),'quantile',probs=c(0.025,0.975))
        }
        K[[k+1]] <- list('beta'=beta,
                         'beta_mean'=b_mean, 'beta_mean_CI'=b_mean_CI,
                         'beta_var'=beta_sd^2, 'beta_var_CI'=beta_sd_CI^2, # provide variance estimates as output
                         'cp'=cp,
                         'cp_mean'=cp_mean, 'cp_mean_CI'=cp_mean_CI,
                         'cp_var'=cp_sd^2,'cp_var_CI'=cp_sd_CI^2) # provide variance estimates as output
      }
      names(K)[k+1] <- paste0('K_', k)
      param_est[[c]] <- list('K_prob'=K_prob,'K'=K)
    }
  }
  names(param_est) <- paste0('Class_', 1:n_class)

  param_est$error_var <- apply(full_out$sigma2_error, c(1), 'mean')

  if(n_cov_outcome_predictive>0){
    outcome_predictive_covariate_alpha <- apply(full_out$outcome_predictive_covariate_alpha, c(1), 'mean')
    param_est$outcome_predictive_covariates <- outcome_predictive_covariate_alpha
    names(param_est$outcome_predictive_covariates) <- outcome_predictive_vars
  }

  if(n_cov_class_predictive>0){
    class_predictive_covariate_lambda <- apply(full_out$class_predictive_covariate_lambda, c(1), 'mean')
    param_est$class_predictive_covariates <- class_predictive_covariate_lambda
    names(param_est$class_predictive_covariates) <- class_predictive_vars

    logistic_intercept <- apply(full_out$logistic_intercept, c(1), 'mean')
    param_est$logistic_intercept <- logistic_intercept
  }

  ## Stop tracking run time
  run_time_total_end <- Sys.time()
  run_time_total <- run_time_total_end - run_time_total_start

  my_results <- list('Convergence'=convergence,
                     'Model_Fit'=model_fit,
                     'Fitted_Values'=y_mean,
                     'Class_Information'=class_info,
                     'Parameter_Estimates'=param_est,
                     'Run_Time'=format(run_time_total))
  if(save_full_chains==TRUE){my_results$Full_MCMC_Chains=full_out}
  if(save_conv_chains==TRUE){my_results$Convergence_MCMC_Chains=mcmc_list}
  class(my_results)='PREM'
  return(my_results)
}
