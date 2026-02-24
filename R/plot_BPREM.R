#' Plot the results of a bivariate piecewise random effects model (BPREM)
#'
#' @description
#' Provides a fitted plot of a BPREM model, as returned by `Bayes_BPREM()`.
#'
#' @param object An object of class "BPREM" (returned by `Bayes_BPREM(...)`).
#' @param xlab X-axis label for the generated plot.
#' @param ylab Y-axis label for the generated plot.
#' @param colors Colors for each class outcome. By default, up to 2 colors are provided in the following order: "blue" (outcome 1), "red" (outcome 2).
#' @param mean_colors Colors for the trajectory defined by the mean parameters for each outcome. By default, up to 5 colors are provided in the following order: "darkblue" (outcome 1), "darkred" (outcome 2).
#' @param legend_pos (optional) Option to change legend position (default = "topright").
#' @param ... (optional) Other parameters to pass to the `plot()` function.
#'
#' @returns No return value.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_bprem)
#' # plot fitted results
#' plot(results_bprem)
#'
#' @import graphics
#'
#' @export
plot.BPREM <- function(object,
                       xlab='X', ylab='Y',
                       colors=NULL, mean_colors=NULL,
                       legend_pos="topright", ...){
  # Setup ----

  data <- object$Data
  id_var <- object$Call$id_var
  time_var <- object$Call$time_var
  y_var <- object$Call$y1_var
  y2_var <- object$Call$y2_var

  if(is.null(colors)) colors <- c('blue','red')
  if(is.null(mean_colors)) mean_colors = c('darkblue','darkred')

  ## outcome data - matrix form
  y <- reshape(data[,c(id_var, time_var, y_var)],
               idvar=id_var,
               timevar=time_var,
               direction='wide')
  y <- unname(as.matrix(y[,names(y)!=id_var]))

  y2 <- reshape(data[,c(id_var, time_var, y2_var)],
                idvar=id_var,
                timevar=time_var,
                direction='wide')
  y2 <- unname(as.matrix(y2[,names(y2)!=id_var]))

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

  # pull fixed effect estimates
  mean_est1 <- object$Parameter_Estimates$Mean[1:4]
  mean_est2 <- object$Parameter_Estimates$Mean[5:8]

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
  legend(legend_pos, lty=1, col=colors[1:2], legend=c(y_var, y2_var))

}
