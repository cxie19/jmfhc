#' Estimate cure rates for subgroups
#' @description Estimates a cure rate for a subgroup of patients with specified characteristics.
#' @param object an object of jmfhc_point_est.
#' @param z_value value(s) of the long-term baseline covariate(s) in the order specified at the argument beta_variable of jmfhc_point_est() function.
#' 
#' @return an estimated cure rate
#' @export
#'
#' @examples
#' result_coef <- jmfhc_point_est(data=jmfhc_dat, event_time="event.time", event_status="event", 
#'                                id="patient.id", beta_variable="trt", gamma_variable="trt",
#'                                fu_measure="measure", fu_time_variable="mes.times")
#' # cure rate for patient subgroup receiving treatment A (trt=0)
#' estcure(object=result_coef,z_value=0)
#' # cure rate for patient subgroup receiving treatment B (trt=1)
#' estcure(object=result_coef,z_value=1)
#' @import cubature
#' @import mvtnorm
#'
estcure <- function(object, z_value){
  coef <- object$coef
  Sigma <- object$re_cov
  length_re <- nrow(Sigma)
  cure_comp <- function(x) {
    X <- matrix(c(x[1],x[2]),nrow=2)
    f <- dmvnorm(x, mean = rep(0,length_re), sigma = Sigma,log=F)
    return(exp(-exp(coef[seq(1+length(z_value)+length_re)]%*%c(1,z_value,X)))*f)
  }
  return(pcubature(cure_comp, lowerLimit = c(-Inf, -Inf), upperLimit = c(+Inf,+Inf))$integral)
}
