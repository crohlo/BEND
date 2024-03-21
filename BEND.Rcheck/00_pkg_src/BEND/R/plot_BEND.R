#' Plot a BEND Model (PREM, CREM, BPREM)
#'
#' @description
#' Generates a "spaghetti plot" of observed longitudinal trajectories for each individual. If the results from a `BEND` function are supplied, the trajectory defined by the mean parameters is shown in bold. If fitting a mixture (`PREMM` or `CI-PREMM`) or bivariate model (`BPREM`), the mean trajectories for classes or outcomes will be distinguished by color.
#'
#' @param data Data frame in long format, where each row describes a measurement occasion for a given individual. It is assumed that each individual has the same number of assigned timepoints (a.k.a., rows).
#' @param id_var Name of column that contains ids for individuals with repeated measures in a longitudinal dataset.
#' @param time_var Name of column that contains the time variable.
#' @param y_var Name of column that contains the outcome variable.
#' @param y2_var (for `BPREM` only) Name of column that contains the second outcome variable.
#' @param results The output of `BEND` model to the data. If results=NULL, only a spaghetti plot of the data will be generated.
#' @param xlab X-axis label for the generated plot.
#' @param ylab Y-axis label for the generated plot.
#' @param colors Colors for each class (`PREMM` or `CI-PREMM`) or outcome (`BPREM`). By default, up to 5 colors are provided in the following order: “blue” (class 1 and outcome 1), “red” (class 2 and outcome 2), “green” (class 3), “gold” (class 4), “gray” (class 5).
#' @param mean_colors Colors for the trajectory defined by the mean parameters for each class (`PREMM` or `CI-PREMM`) or outcome (`BPREM`). By default, up to 5 colors are provided in the following order: “darkblue” (class 1 and outcome 1), “darkred” (class 2 and outcome 2), “darkgreen” (class 3), “gold4” (class 4), “darkgray” (class 5).
#' @param legend_pos (optional) Option to change legend position (default = "topright").
#' @param ... (optional) Other parameters to pass to the `plot()` function.
#'
#' @returns No return value, called to generate plot.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load simulated data
#' data(SimData_PREM)
#' # plot observed data
#' plot_BEND(data = SimData_PREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y")
#' # load fitted model results
#' data(results_prem)
#' # plot fitted results
#' plot_BEND(data = SimData_PREM,
#'           id_var = "id",
#'           time_var = "time",
#'           y_var = "y",
#'           results = results_prem)
#'
#' @import graphics
#'
#' @export
plot_BEND <- function(data,
                      id_var, time_var, y_var,
                      y2_var=NULL,
                      results=NULL,
                      xlab='X', ylab='Y',
                      colors=NULL, mean_colors=NULL,
                      legend_pos="topright", ...){
  # Setup ----

  if(is.null(colors)) colors <- c('blue','red','green','gold','gray')
  if(is.null(mean_colors)) mean_colors = c('darkblue','darkred','darkgreen','gold4','darkgray')

  ## outcome data - matrix form
  y <- reshape(data[,c(id_var, time_var, y_var)],
               idvar=id_var,
               timevar=time_var,
               direction='wide')
  y <- unname(as.matrix(y[,names(y)!=id_var]))

  # For BPREM only
  if(!is.null(y2_var)){
    y2 <- reshape(data[,c(id_var, time_var, y2_var)],
                  idvar=id_var,
                  timevar=time_var,
                  direction='wide')
    y2 <- unname(as.matrix(y2[,names(y2)!=id_var]))
  }

  ## time data - matrix form
  ## should be the same dimensions as y
  x <- matrix(data[,c(time_var)],
              byrow=TRUE,
              nrow=dim(y)[1],
              ncol=dim(y)[2])

  ## Define relevant variables
  n_subj <- nrow(y)
  n_time <- ncol(y)
  xvec <- seq(min(x), max(x), length.out=100)

  # Observed Plot -----
  if(is.null(results)){
    plot(x[1, !is.na(y[1,])],
         y[1, !is.na(y[1,])],
         type = "l",
         ylim = c(min(y,na.rm=TRUE), max(y,na.rm=TRUE)),
         xlim= c(min(x,na.rm=TRUE), max(x,na.rm=TRUE)),
         xlab = xlab,
         ylab = ylab, ...)
    for(i in 2:n_subj){
      lines(x[i, !is.na(y[i,])],
            y[i, !is.na(y[i,])],
            type = "l")
    }

    # BPREM only
    if(!is.null(y2_var)){
      for(i in 2:n_subj){
        lines(x[i, !is.na(y[i,])],
              y[i, !is.na(y[i,])],
              type = "l",
              col=colors[1])
      }
      for(i in 1:n_subj){
        lines(x[i, !is.na(y2[i,])],
              y2[i, !is.na(y2[i,])],
              type = "l",
              col = colors[2])
      }
      legend(legend_pos, lty=1, col=colors[1:2], legend=c("Outcome 1", "Outcome 2"))
    }

    return('Observed trajectories')
  }

  # Fitted Plot -----
  if(!is.null(results)){

    ## PREM -----
    if(inherits(results, "PREM")){

      # determine number of classes
      n_class <- length(unique(results$Class_Information$class_membership))

      # determine number of changepoints in each class (based on final model results)
      changepoints <- c()
      class_data <- data.frame()
      for(i in 1:n_class){
        class_num <- paste0("Class_", i)
        changepoints[i] <- which.max(results$Parameter_Estimates[[class_num]]$K_prob)-1
      }
      max_cp <- max(changepoints)

      # determine who is in each class
      class_list <- list()
      for(i in 1:n_class){
        class_list[[i]] <- c(1:n_subj)[results$Class_Information$class_membership==i]
      }

      # class mean estimates
      class_means <- list()
      for(i in 1:n_class){
        class_num <- paste0("Class_", i)
        cp_num <- paste0("K_", changepoints)[i]
        k <- changepoints[i]

        if(k==0) I <- rep(0, max_cp)
        if(k>0) I <- c(rep(1,k), rep(0,max_cp-k+2))
        # i = k+1
        int <- results$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[1]
        slope1 <- results$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[2]

        slope_cp <- rep(0,max_cp)
        cp <- rep(0,max_cp)
        if(k>0){
          slope_cp[1:k] <- results$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[3:(k+2)]
          cp[1:k] <- results$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_mean[1:k]
        }

        class_means[[i]] <- rep(0,100)
        for(j in 1:100){
          temp <- int + slope1*xvec[j]
          if(k>0){
            for(l in 1:k){
              temp <- temp + slope_cp[l]*(max(0, xvec[j]-cp[l]))}
          }
          class_means[[i]][j] <- temp
        }
      }

      plot(NULL, NULL,
           xlim = c(min(x), max(x)),
           ylim = c(min(y,na.rm=TRUE), max(y,na.rm=TRUE)),
           ylab = ylab,
           xlab = xlab, ...)
      for(i in 1:n_class){
        for(j in class_list[[i]]){
          points(x[j,!is.na(y[j,])],
                 y[j,!is.na(y[j,])],
                 type = "l",
                 col=colors[i])
        }
        points(xvec,
               class_means[[i]],
               type='l',
               col=mean_colors[i],
               lwd=4)
      }
      legend(legend_pos, lty=1, col=colors[1:n_class], legend=paste0("Class ", 1:n_class))
    }

    ## CREM -----
    if(inherits(results, "CREM")){

      # determine form
      form <- results$Functional_Form

      # determine number of parameters
      if(form=="linear")                           n_param <- 2
      if(form=="quadratic" | form=="exponential")  n_param <- 3
      if(form=="piecewise")                        n_param <- 4

      # pull fixed effect estimates
      mean_est <- results$Parameter_Estimates$Mean[1:n_param]

      # define functional form equation

      fit_form_eq <- function(x,est){
        if(form=="linear")       return(est[1] + est[2]*x)
        if(form=="quadratic")    return(est[1] + est[2]*x + est[3]*(x^2))
        if(form=="exponential")  return(est[1] + est[2]*(1-exp(-est[3]*x)))
        if(form=="piecewise")    return(est[1] + est[2]*x + est[3]*(max(0, x-est[4])))
      }

      mean_traj <- rep(0,100)
      for(i in 1:100){
        mean_traj[i] <- fit_form_eq(xvec[i], mean_est)
      }

      plot(x[1, !is.na(y[1,])],
           y[1, !is.na(y[1,])],
           type = "l",
           col = "grey",
           ylim = c(min(y,na.rm=TRUE), max(y,na.rm=TRUE)),
           xlim= c(min(x,na.rm=TRUE), max(x,na.rm=TRUE)),
           xlab = xlab,
           ylab = ylab, ...)
      for(i in 2:n_subj){
        lines(x[i, !is.na(y[i,])],
              y[i, !is.na(y[i,])],
              type = "l",
              col = "grey")
      }
      points(xvec,
             mean_traj,
             type='l',
             col="black",
             lwd=4)
    }

    ## BPREM -----
    if(inherits(results, "BPREM")){

      # pull fixed effect estimates
      mean_est1 <- results$Parameter_Estimates$Mean[1:4]
      mean_est2 <- results$Parameter_Estimates$Mean[5:8]

      # define piecewise function
      pw_eq <- function(x,est){
        return(est[1] + est[2]*x + est[3]*(max(0, x-est[4])))
      }

      mean_traj1 <- rep(0,100)
      mean_traj2 <- rep(0,100)
      for(i in 1:100){
        mean_traj1[i] <- pw_eq(xvec[i], mean_est1)
        mean_traj2[i] <- pw_eq(xvec[i], mean_est2)
      }

      plot(x[1, !is.na(y[1,])],
           y[1, !is.na(y[1,])],
           type = "l",
           col = "grey",
           ylim = c(min(y,na.rm=TRUE), max(y,na.rm=TRUE)),
           xlim= c(min(x,na.rm=TRUE), max(x,na.rm=TRUE)),
           xlab = xlab,
           ylab = ylab, ...)
      for(i in 2:n_subj){
        lines(x[i, !is.na(y[i,])],
              y[i, !is.na(y[i,])],
              type = "l",
              col = colors[1])
      }
      for(i in 1:n_subj){
        lines(x[i, !is.na(y2[i,])],
              y2[i, !is.na(y2[i,])],
              type = "l",
              col = colors[2])
      }
      points(xvec,
             mean_traj1,
             type='l',
             col=mean_colors[1],
             lwd=4)
      points(xvec,
             mean_traj2,
             type='l',
             col=mean_colors[2],
             lwd=4)
      legend(legend_pos, lty=1, col=colors[1:2], legend=c("Outcome 1", "Outcome 2"))
    }

  return('Observed trajectories with fitted results')

  }
}
