% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jmfhc_se_est.R
\name{jmfhc_se_est}
\alias{jmfhc_se_est}
\title{Obtain the Standard Error Estimation of the Joint Cure Model with Flexible Hazards Ratios Over Time}
\usage{
jmfhc_se_est(
  data,
  event_time,
  event_status,
  id,
  beta_variable,
  gamma_variable = NULL,
  fu_measure,
  fu_time_fixed_variable,
  fu_time_random_variable,
  baseline_var_lmm = NULL,
  max.int = 200,
  no_cores = 5,
  no_cores_jackknife = 5,
  absconverge.par = 0.001,
  relconverge.par = 0.002,
  absconverge.F0t = 0.002,
  relconverge.F0t = 0.005
)
}
\arguments{
\item{data}{a data.frame in a long format with no missing values contains patients' id, observed survival time, event status, any continuous variable(s) and/or dummy variable(s) of
categorical variable(s) at baseline, and longitudinal measurements with the corresponding measurement time points. The baseline categorical variables with more than two levels need to be
converted into dummy variables before fitting the function. The longitudinal measurements could be in the original form or the transformed form(s) for the longitudinal submodel fitting.
The data set needs to contain the original form of measurement time points in the same unit of observed survival times. If transformed forms of measurement time points are used for
longitudinal submodel fitting, the data set also needs to have columns for the transformed forms of measurement time points. All patients need to have at least two measurements
including the baseline.}

\item{event_time}{the variable name of observed survival time in the data.}

\item{event_status}{the variable name of event status in the data.}

\item{id}{the variable name of patient id in the data.}

\item{beta_variable}{the variable name(s) defined as long-term covariate(s) in the model, which cannot be NULL.}

\item{gamma_variable}{the variable name(s) defined as short-term covariate(s) in the model. By default gamma_variable = NULL.
If there is no defined short-term covariate, the cure submodel of the joint model becomes a proportional hazards cure submodel.}

\item{fu_measure}{the variable name of repeatedly measured outcome (e.g., a biomarker) in the original form or any transformation
used as the outcome of the longitudinal submodel (linear mixed-effects model).}

\item{fu_time_fixed_variable}{the variable name(s) of functions of measurement times used as the fixed effect(s) in
the longitudinal submodel.}

\item{fu_time_random_variable}{the variable name(s) of functions of measurement times used as the random effect(s) in
the longitudinal submodel.}

\item{baseline_var_lmm}{the baseline variable name(s) in the longitudinal submodel. By default baseline_var_lmm = NULL.}

\item{max.int}{maximum number of iterations. The default is 200.}

\item{no_cores}{the number of cores used during each point estimation for one resampled data set. The default is 5.}

\item{no_cores_jackknife}{the number of cores used to estimate jackknife resampled data sets. The default is 5.}

\item{absconverge.par}{absolute difference \eqn{(current-previous)} for regression parameters as the convergence criteria. The default is 1e-3.}

\item{relconverge.par}{relative difference \eqn{[(current-previous)/previous]} for regression parameters as the convergence criteria. The default is 2e-3.}

\item{absconverge.F0t}{the average of absolute differences for \eqn{F_0(t)} as the convergence criteria. The default is 2e-3.}

\item{relconverge.F0t}{the average of relative differences for \eqn{F_0(t)} as the convergence criteria. The default is 5e-3.}
}
\value{
a list containing results of regression parameters' standard errors and point estimation results from jackknife resampled datasets.
The following items std.err and result.jackknife are returned.
\item{std.err}{standard errors of regression parameters (beta and gamma)}
\item{result.jackknife}{a matrix containing point estimates of regression parameters from joint models for jackknife resampled datasets}
}
\description{
Fits the joint cure model and estimates standard errors of regression parameters via a jackknife resampling method
by removing one patient's record each time.
Tied failure times are dealt with Breslow approximation.
}
\examples{
data(jmfhc_dat)
result_se <- jmfhc_se_est(data=jmfhc_dat, event_time="event.time", event_status="event",
                          id="patient.id", beta_variable="trt", gamma_variable="trt",
                          fu_measure="measure", fu_time_fxied_variable="mes.times",
                          fu_time_random_variable="mes.times")
result_se <- jmfhc_se_est(data=jmfhc_dat, event_time="event.time", event_status="event",
                          id="patient.id", beta_variable="trt", fu_measure="measure",
                          fu_time_fxied_variable="mes.times",fu_time_random_variable="mes.times")
result$coef
}
