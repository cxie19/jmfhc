#' A toy data set used for modeling JMFHC
#'
#' This toy data set with a cure fraction contains 500 patients' time-to-event outcomes, basline information,
#' and a biomarker's repeated measured values.
#'
#'@format A data frame with 500 observations on the following 6 variables.
#' \describe{
#'     \item{patient.id}{patient id}
#'     \item{mes.time}{measurement time points for the biomarker}
#'     \item{measure}{repeatedly measure biomarker values}
#'     \item{trt}{a binary variable with 0 for treatment A and 1 for treatment B}
#'     \item{event.time}{observed survival time}
#'     \item{event}{an event indicator with 1 for dead and 0 for censored}
#'     }
#'
#'@source{Generated from JMFHC to serve as an example.}
#'
#'@examples
#' data(jmfhc_dat)
"jmfhc_dat"
