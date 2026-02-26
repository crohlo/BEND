#' Plot the results of a crossed random effects model (CREM)
#'
#' @description
#' Provides a fitted plot of a CREM model, as returned by `Bayes_CREM()`.
#'
#' @param x An object of class "CREM" (returned by `Bayes_CREM(...)`).
#' @param xlab X-axis label for the generated plot.
#' @param ylab Y-axis label for the generated plot.
#' @param colors Color for observed trajectories (optional). Default is "grey".
#' @param mean_colors Colors for the trajectory defined by the mean parameters for each outcome (optional). Default is "black".
#' @param legend_pos (optional) Option to change legend position (default = "topright").
#' @param ... (optional) Other parameters to pass to the `plot()` function.
#'
#' @returns No return value.
#'
#' @author Corissa T. Rohloff
#'
#' @examples
#' # load fitted model results
#' data(results_pcrem)
#' # plot fitted results
#' plot(results_pcrem)
#'
#' @import graphics
#'
#' @export
plot.CREM <- function(x,
                      xlab='X', ylab='Y',
                      colors=NULL, mean_colors=NULL,
                      legend_pos="topright", ...){
  # Setup ----

  if(is.null(colors)) colors <- c('grey')
  if(is.null(mean_colors)) mean_colors = c('black')

  object <- x

  data <- object$Data
  id_var <- object$Call$ind_id_var
  time_var <- object$Call$time_var
  y_var <- object$Call$y_var

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

  # determine form
  form <- object$Functional_Form

  # determine number of parameters
  if(form=="linear")                           n_param <- 2
  if(form=="quadratic" | form=="exponential")  n_param <- 3
  if(form=="piecewise")                        n_param <- 4

  # pull fixed effect estimates
  mean_est <- object$Parameter_Estimates$Mean[1:n_param]

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
       col = colors,
       ylim = c(min(y,na.rm=TRUE), max(y,na.rm=TRUE)),
       xlim= c(min(x,na.rm=TRUE), max(x,na.rm=TRUE)),
       xlab = xlab,
       ylab = ylab, ...)
  for(i in 2:n_subj){
    lines(x[i, !is.na(y[i,])],
          y[i, !is.na(y[i,])],
          type = "l",
          col = colors)
  }
  points(xvec,
         mean_traj,
         type='l',
         col=mean_colors,
         lwd=4)

}
