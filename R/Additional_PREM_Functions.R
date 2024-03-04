# Additional Functions ----

## Permute Changepoint Labels ----
Permute_CP <- function(output, n_class, max_cp, max_beta){
  # output = output from JAGS
  # n_class = number of latent clases
  # max_cp = maximum number of changepoints
  # max_beta = maximum number of intercept + slope parameters

  ### Setup -----
  n_iter = dim(output$sigma2_error)[2]
  n_chain = dim(output$sigma2_error)[3]

  ### Permute changepoint labels, if necessary, so they are ordered -----
  for(i in 1:n_iter){
    for(j in 1:n_chain){
      for(k in 1:n_class){
        if(output$n_cp[k,i,j]>0){
          ind <- c(1:max_cp)[as.logical(output$cp_indicator[k,,i,j])]
          my_order <- order(output$cp[k,ind,i,j])
          output$cp_indicator[k,,i,j] <- c(output$cp_indicator[k,ind[my_order],i,j],output$cp_indicator[k,-ind,i,j])
          output$cp[k,,i,j] <- c(output$cp[k,ind[my_order],i,j],output$cp[k,-ind,i,j])
          output$beta[k,3:max_beta,i,j] <- c(output$beta[k,ind[my_order]+2,i,j],output$beta[k,c(3:max_beta)[-ind],i,j])
        }
      }
    }
  }
  return(output)
}


## Get Mode (for initializing number of changepoints) ----
getmode <- function(v){
  # v = vector

  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


## Realign ECR Function ----
Realign_ECR <- function(output, n_class, model=NULL){
  # output = output from JAGS
  # n_class = number of latent clases
  # model = model type ("PREMM", "CI_PREMM_Full", "CI_PREMM_Class_Predictive" or "CI_PREMM_Outcome_Predictive")

  ### Warnings -----
  if(is.null(model)) stop('Specify model type. Must be \'PREMM\', \'CI_PREMM_Full\', \'CI_PREMM_Class_Predictive\', or \'CI_PREMM_Outcome_Predictive\'')

  ### Setup -----
  n_subj <- dim(output$class)[1]
  n_iter <- dim(output$class)[2]
  n_chain <- dim(output$class)[3]
  if(model=="CI_PREMM_Full" | model=="CI_PREMM_Class_Predictive"){class_pred <- TRUE} else {class_pred <- FALSE}
  ## For class_pred == TRUE
  if(class_pred == TRUE){output$class_orig <- output$class}

  ### Align for the first time -----
  class_matrix <- matrix(nrow=n_iter*n_chain, ncol=n_subj)
  for(i in 1:n_subj){
    class_matrix[,i] <- as.vector(output$class[i,,])
  }
  permutes <- label.switching::ecr.iterative.1(class_matrix, n_class)$permutations

  perm_chain <- list()
  perm_chain[[1]] <- permutes[1:n_iter,]
  perm_chain[[2]] <- permutes[(n_iter+1):(2*n_iter),]
  perm_chain[[3]] <- permutes[(2*n_iter+1):(3*n_iter),]
  for(i in 1:n_chain){
    for(j in 1:n_iter){
      # this loop will ignore parameters that don't exist in output
      ## Initialization-Only
      output$beta[,,j,i] <- output$beta[perm_chain[[i]][j,],,j,i]
      output$cp[,,j,i] <- output$cp[perm_chain[[i]][j,],,j,i]
      ## PREMM only
      output$class_prob[,j,i] <- output$class_prob[perm_chain[[i]][j,],j,i]
      ## Full Model
      output$beta_mean[,,j,i] <- output$beta_mean[perm_chain[[i]][j,],,j,i]
      output$beta_sd[,,j,i] <- output$beta_sd[perm_chain[[i]][j,],,j,i]
      output$cp_mean[,,j,i] <- output$cp_mean[perm_chain[[i]][j,],,j,i]
      output$cp_sd[,,j,i] <- output$cp_sd[perm_chain[[i]][j,],,j,i]
      ## ALL
      output$cp_indicator[,,j,i] <- output$cp_indicator[perm_chain[[i]][j,],,j,i]
      output$n_cp[,j,i] <- output$n_cp[perm_chain[[i]][j,],j,i]
      inv_perm <- order(perm_chain[[i]][j,])
      output$class[,j,i] <- inv_perm[output$class[,j,i]]
    }
  }

  ### Align for a second time (for better accuracy) -----
  class_matrix <- matrix(nrow=n_iter*n_chain, ncol=n_subj)
  for(i in 1:n_subj){
    class_matrix[,i] <- as.vector(output$class[i,,])
  }
  permutes <- label.switching::ecr.iterative.1(class_matrix, n_class)$permutations

  perm_chain <- list()
  perm_chain[[1]] <- permutes[1:n_iter,]
  perm_chain[[2]] <- permutes[(n_iter+1):(2*n_iter),]
  perm_chain[[3]] <- permutes[(2*n_iter+1):(3*n_iter),]
  for(i in 1:n_chain){
    for(j in 1:n_iter){
      # this loop will ignore parameters that don't exist in output
      ## Initialization-Only
      output$beta[,,j,i] <- output$beta[perm_chain[[i]][j,],,j,i]
      output$cp[,,j,i] <- output$cp[perm_chain[[i]][j,],,j,i]
      ## PREMM only
      output$class_prob[,j,i] <- output$class_prob[perm_chain[[i]][j,],j,i]
      ## Full Model
      output$beta_mean[,,j,i] <- output$beta_mean[perm_chain[[i]][j,],,j,i]
      output$beta_sd[,,j,i] <- output$beta_sd[perm_chain[[i]][j,],,j,i]
      output$cp_mean[,,j,i] <- output$cp_mean[perm_chain[[i]][j,],,j,i]
      output$cp_sd[,,j,i] <- output$cp_sd[perm_chain[[i]][j,],,j,i]
      ## ALL
      output$cp_indicator[,,j,i] <- output$cp_indicator[perm_chain[[i]][j,],,j,i]
      output$n_cp[,j,i] <- output$n_cp[perm_chain[[i]][j,],j,i]
      inv_perm <- order(perm_chain[[i]][j,])
      output$class[,j,i] <- inv_perm[output$class[,j,i]]
    }
  }

  ### Invert logistic regression models if the class was switched -----
  if(class_pred == TRUE){
    n_cov <- length(output$class_predictive_covariate_lambda[,1,1])
    for(i in 1:n_iter){
      for(j in 1:n_chain){
        for(k in 1:n_cov){

          if(output$class[1,i,j] != output$class_orig[1,i,j]){
            cp_class1 = output$cp_indicator[1,,i,j]
            cp_class2 = output$cp_indicator[2,,i,j]

            output$cp_indicator[1,,i,j] = cp_class2
            output$cp_indicator[2,,i,j] = cp_class1

            output$logistic_intercept[1,i,j] = -output$logistic_intercept[1,i,j]
            output$class_predictive_covariate_lambda[k,i,j] = -output$class_predictive_covariate_lambda[k,i,j]
          }

        }
      }
    }
  }
  return(output)
}

