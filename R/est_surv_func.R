#' Estimate survival functions for subgroups
#' @description Estimates the survival function for a subgroup of patients with specified characteristics.
#' @param object an object of jmfhc_point_est.
#' @param z_value value(s) of the long-term baseline covariate(s) in the order specified at the argument beta_variable of jmfhc_point_est() function.
#' @param x_value value(s) of the short-term baseline covariate(s) in the order specified at the argument gamma_variable of jmfhc_point_est() function.
#' If gamma is set as NULL in the jmfhc_point_est() function, then there is no need to put values here. The default is NULL.
#' @param event_time the variable name of observed survival time in the data.
#' @param n.random the number of randomly generated random effects based on the estimated covariance matrix for random effects. The default is 1000.
#'
#' @return a data.frame with sorted survival times and the corresponding estimated survival probabilities.
#' @export
#'
#' @examples
#' result_coef <- jmfhc_point_est(data=jmfhc_dat, event_time="event.time", event_status="event", 
#'                                id="patient.id", beta_variable="trt", gamma_variable="trt",
#'                                fu_measure="measure", fu_time_variable="mes.times")
#' # estimated survival function for patient subgroup receiving treatment A
#' survival_trt0 <- est_surv_func(object=result_coef,z_value=0,x_value=0,event_time="event.time")
#' # estimated survival function for patient subgroup receiving treatment B
#' survival_trt1 <- est_surv_func(object=result_coef,z_value=1,x_value=1,event_time="event.time")
#' @import mvtnorm

est_surv_func <- function(object,z_value,x_value=NULL,event_time,n.random=1000){
  
  set.seed(1)
  
  coef <- object$coef
  Sigma <- object$re_cov
  length_re <- nrow(Sigma)
  betas <- coef[seq(1+length(z_value)+length_re)]
  dat_base <- object$dat_baseline
  F0t <- dat_base$base_cdf[order(dat_base[,event_time])]
  mvn <- mvrnorm(n.random, mu = rep(0,length_re), Sigma = Sigma)
  
  if (!is.null(x_value)){
    gammas <- coef[(1+length(z_value)+length_re)+seq(x_value)]
    est_survival <- sapply(seq(n.random),function(x)
      exp(-rep(exp(betas%*%c(1,z_value,mvn[x,])),length(F0t))*(F0t^rep(exp(gammas%*%x_value),length(F0t)))))
  }else{
    est_survival <- sapply(seq(n.random),function(x)
      exp(-rep(exp(betas%*%c(1,z_value,mvn[x,])),length(F0t))*F0t))
  }
  
  survival_mean <- apply(est_survival,1,mean)
  
  return(data.frame(time=sort(dat_base[,event_time]),survival=survival_mean))
  
}

