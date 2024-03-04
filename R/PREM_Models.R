# Bayes_PREM JAGS Models

# PREM (One class Models) ----

## Binomial CP ----

### Fixed ----

model_prem_binomial_fixed <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
      }
    }
  }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

### Scale Unif ----

model_prem_binomial_scaleunif <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Scale HC ----

model_prem_binomial_scalehc <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

## Uniform CP ----

### Fixed ----

model_prem_uniform_fixed <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
       x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

### Scale Unif ----

model_prem_uniform_scaleunif <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Scale HC ----

model_prem_uniform_scalehc <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

## Binomial CP (Covariates) ----

### Fixed ----

model_prem_cov_binomial_fixed <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
       x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
      }
    }
  }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

### Scale Unif ----

model_prem_cov_binomial_scaleunif <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Scale HC ----

model_prem_cov_binomial_scalehc <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"


## Uniform CP (Covariates) ----

### Fixed ----

model_prem_cov_uniform_fixed <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
       x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

### Scale Unif ----

model_prem_cov_uniform_scaleunif <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Scale HC ----

model_prem_cov_uniform_scalehc <- "model{
  for(i in 1:n_subj){
  class[i] <- 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

# PREMM (No Covariates) ----

## Binomial CP ----

### Fixed ----

model_premm_binomial_fixed <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
       x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
      }
    }
  }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

### Scale Unif ----

model_premm_binomial_scaleunif <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Scale HC ----

model_premm_binomial_scalehc <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

## Uniform CP ----

### Fixed ----

model_premm_uniform_fixed <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
       x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
 }
"

### Scale Unif ----

model_premm_uniform_scaleunif <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Scale HC ----

model_premm_uniform_scalehc <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

# CI-PREMM ----

## Full ----

### Binomial CP ----

#### Fixed ----

model_cipremm_full_binomial_fixed <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
      }
    }
  }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

#### Scale Unif ----

model_cipremm_full_binomial_scaleunif <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

#### Scale HC ----

model_cipremm_full_binomial_scalehc <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Uniform CP ----

#### Fixed ----

model_cipremm_full_uniform_fixed <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

#### Scale Unif ----

model_cipremm_full_uniform_scaleunif <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

#### Scale HC ----

model_cipremm_full_uniform_scalehc <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

## Class Predictive - Only (CPO) ----

### Binomial CP ----

#### Fixed ----

model_cipremm_cpo_binomial_fixed <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
      }
    }
  }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

#### Scale Unif ----

model_cipremm_cpo_binomial_scaleunif <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

#### Scale HC ----

model_cipremm_cpo_binomial_scalehc <- "model{
  for(i in 1:n_subj){
    logit(logistic_class_prob[i]) <- logistic_intercept +
      sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
    class_bern[i] ~ dbern(logistic_class_prob[i])
    class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
      y[i,j] ~ dnorm(mu_y[i,j], tau_y)
      mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                   sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
      beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
      cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Uniform CP ----

#### Fixed ----

model_cipremm_cpo_uniform_fixed <- "model{
  for(i in 1:n_subj){
  logit(logistic_class_prob[i]) <- logistic_intercept +
                              sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
  class_bern[i] ~ dbern(logistic_class_prob[i])
  class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

#### Scale Unif ----

model_cipremm_cpo_uniform_scaleunif <- "model{
  for(i in 1:n_subj){
    logit(logistic_class_prob[i]) <- logistic_intercept +
      sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
    class_bern[i] ~ dbern(logistic_class_prob[i])
    class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
      y[i,j] ~ dnorm(mu_y[i,j], tau_y)
      mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                   sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
      beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
      cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

#### Scale HC ----

model_cipremm_cpo_uniform_scalehc <- "model{
  for(i in 1:n_subj){
    logit(logistic_class_prob[i]) <- logistic_intercept +
      sum(class_predictive_covariate_lambda[1:n_cov_class_predictive]*class_predictive_covariates[i,1:n_cov_class_predictive])
    class_bern[i] ~ dbern(logistic_class_prob[i])
    class[i] <- class_bern[i] + 1
    for(j in 1:n_time){
      y[i,j] ~ dnorm(mu_y[i,j], tau_y)
      mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
        sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
      beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
      cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the logistic intercept
  logistic_intercept ~ dnorm(0, 0.1)

  # prior distribution of the class predictive covariates
  for(p in 1:n_cov_class_predictive){
  class_predictive_covariate_lambda[p] ~ dnorm(0, 0.1)
  }

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

## Outcome Predictive - Only (OPO) ----

### Binomial CP ----

#### Fixed ----

model_cipremm_opo_binomial_fixed <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
      }
    }
  }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

#### Scale Unif ----

model_cipremm_opo_binomial_scaleunif <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
   for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

#### Scale HC ----

model_cipremm_opo_binomial_scalehc <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
   for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    cp_indicator[c,k] ~ dbern(binom_prob)

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

### Uniform CP ----

#### Fixed ----

model_cipremm_opo_uniform_fixed <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
    for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[class[i],1] + beta[class[i],2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[class[i],3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])
    # prior distribution of fixed effect of changepont
    cp[c,k] ~ dunif(min_cp_mean, max_cp_mean)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[c,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])

  # prior distribution of the fixed parameters
  beta[c,1:max_beta] ~ dmnorm(mean, prec)
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y
  }
"

#### Scale Unif ----

model_cipremm_opo_uniform_scaleunif <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
   for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dunif(0, (max_time-min_time)/4)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dunif(0, (2/prec[q,q])^0.5)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

#### Scale HC ----

model_cipremm_opo_uniform_scalehc <- "model{
  for(i in 1:n_subj){
  class[i] ~ dcat(class_prob)
   for(j in 1:n_time){
    y[i,j] ~ dnorm(mu_y[i,j], tau_y)
    mu_y[i,j] <- beta[i,1] + beta[i,2]*x[i,j] +
                 sum(cp_indicator[class[i],1:max_cp]*beta[i,3:max_beta]*x_cp[class[i],1:max_cp,i,j]) +
                 inprod(outcome_predictive_covariate_alpha[], outcome_predictive_covariates[i,j,])
    }
  }

  # prior distribution of the random coefficients
  for(i in 1:n_subj){
    for(q in 1:max_beta){
    beta[i,q] ~ dnorm(beta_mean[class[i],q], beta_tau[class[i],q])
    }
    for(k in 1:max_cp){
    cp[i,k] ~ dnorm(cp_mean[class[i],k], cp_prec[class[i],k]) T(min_time,max_time)
    }
  }

  # prior distribution of the changepoint
  for(c in 1:n_class){
    for(k in 1:max_cp){
    # prior distribution of number of changepoints
    Temp[c,k] ~ dbern(aux_prob[k])
    cp_indicator[c,k] <- prod(Temp[c,1:k])

    cp_mean[c,k] ~ dunif(min_cp_mean, max_cp_mean)
    cp_prec[c,k] <- 1/(cp_sd[c,k]^2)
    # prior distribution of the changepoint variance
    cp_sd[c,k] ~ dt(0, 1/(((max_time-min_time)/4)/tan(0.45*3.1416))^2, 1) T(0,)
      for(j in 1:n_time){
        for(i in 1:n_subj){
        x_cp[c,k,i,j] <- max(0, x[i,j]-cp[i,k])
        }
      }
    }

  n_cp[c] <- sum(cp_indicator[c,1:max_cp])
  beta_mean[c,1:max_beta] ~ dmnorm(mean, prec)

  # prior distribution of random-effects variances
    for(q in 1:max_beta){
    beta_sd[c,q] ~ dt(0, 1/(((2/prec[q,q])^0.5)/tan(0.45*3.1416))^2, 1) T(0,)
    beta_tau[c,q] <- 1/(beta_sd[c,q]^2)
    }
  }

  # prior distribution of the outcome predictive covariates
  for(l in 1:n_cov_outcome_predictive){
  outcome_predictive_covariate_alpha[l] ~ dnorm(0, 0.001)
  }

  # prior distribution of the class probability
  class_prob ~ ddirch(alpha)

  # prior distribution of the precision for y
  tau_y ~ dgamma(0.001, 0.001)
  sigma2_error <- 1/tau_y

  for_conv[1:n_params] <- c(rep(0, n_params))
  }
"

# CI-PREMM Model Lists ----
## Full ----
cipremm_mods_full = list(model_cipremm_full_binomial_scaleunif,
                         model_cipremm_full_binomial_scalehc,
                         model_cipremm_full_uniform_scaleunif,
                         model_cipremm_full_uniform_scalehc)

## CPO ----
cipremm_mods_cpo = list(model_cipremm_cpo_binomial_scaleunif,
                        model_cipremm_cpo_binomial_scalehc,
                        model_cipremm_cpo_uniform_scaleunif,
                        model_cipremm_cpo_uniform_scalehc)

## OPO ----
cipremm_mods_opo = list(model_cipremm_opo_binomial_scaleunif,
                        model_cipremm_opo_binomial_scalehc,
                        model_cipremm_opo_uniform_scaleunif,
                        model_cipremm_opo_uniform_scalehc)

