#' Plot the results of a piecewise random effects model (PREM)
#'
#' @description
#' Provides a fitted plot of a PREM model, as returned by `Bayes_PREM()`.
#'
#' @param object An object of class "PREM" (returned by `Bayes_PREM(...)`).
#' @param xlab X-axis label for the generated plot.
#' @param ylab Y-axis label for the generated plot.
#' @param colors Colors for each class (`PREMM` or `CI-PREMM`). By default, up to 5 colors are provided in the following order: "blue" (class 1), "red" (class 2), "green" (class 3), "gold" (class 4), "gray" (class 5).
#' @param mean_colors Colors for the trajectory defined by the mean parameters for each class (`PREMM` or `CI-PREMM`). By default, up to 5 colors are provided in the following order: "darkblue" (class 1), "darkred" (class 2), "darkgreen" (class 3), "gold4" (class 4), "darkgray" (class 5).
#' @param legend_pos (optional) Option to change legend position (default = "topright").
#' @param ... (optional) Other parameters to pass to the `plot()` function.
#'
#' @returns No return value.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_prem)
#' # plot fitted results
#' plot(results_prem)
#'
#' @import graphics
#'
#' @export
plot.PREM <- function(object,
                      xlab='X', ylab='Y',
                      colors=NULL, mean_colors=NULL,
                      legend_pos="topright", ...){
  # Setup ----

  data <- object$Data
  id_var <- object$Call$id_var
  time_var <- object$Call$time_var
  y_var <- object$Call$y_var

  if(is.null(colors)) colors <- c('blue','red','green','gold','gray')
  if(is.null(mean_colors)) mean_colors = c('darkblue','darkred','darkgreen','gold4','darkgray')

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

  ## Define relevant variables
  n_subj <- nrow(y)
  n_time <- ncol(y)
  xvec <- seq(min(x), max(x), length.out=100)

  # Fitted Plot -----


  # determine number of classes
  n_class <- length(unique(object$Class_Information$class_membership))

  # determine number of changepoints in each class (based on final model results)
  changepoints <- c()
  class_data <- data.frame()
  for(i in 1:n_class){
    class_num <- paste0("Class_", i)
    changepoints[i] <- which.max(object$Parameter_Estimates[[class_num]]$K_prob)-1
  }
  max_cp <- max(changepoints)

  # determine who is in each class
  class_list <- list()
  for(i in 1:n_class){
    class_list[[i]] <- c(1:n_subj)[object$Class_Information$class_membership==i]
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
    int <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[1]
    slope1 <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[2]

    slope_cp <- rep(0,max_cp)
    cp <- rep(0,max_cp)
    if(k>0){
      slope_cp[1:k] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$beta_mean[3:(k+2)]
      cp[1:k] <- object$Parameter_Estimates[[class_num]]$K[[cp_num]]$cp_mean[1:k]
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
