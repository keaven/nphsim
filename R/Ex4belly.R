#' Time-to-event data example 4 for non-proportional hazards working group
#'
#' Survival objects reverse-engineered datasets from published Kaplan-Meier
#' curves. 
#' Individual trials are de-identified since the data are only
#' approximations of the actual data.
#' Data are intended to evaluate methods and designs for trials where
#' non-proportional hazards may be anticipated for outcome data.  
#'
#' @docType data
#'
#' @usage data(Ex4belly)
#'
#' @format Data frame with 4 variables:
#' \describe{ 
#' \item{id}{sequential numbering of unique identifiers}
#' \item{month}{time-to-event}
#' \item{event}{1 for event, 0 for censored} 
#' \item{trt}{1 for experimental, 0 for control}
#' }
#' 
#' @keywords datasets
#'
#' @references TBD
#' 
#' @seealso \code{\link{Ex1delayedEffect}}, \code{\link{Ex2delayedEffect}}, \code{\link{Ex3belly}},
#' \code{\link{Ex5widening}}, \code{\link{Ex6crossing}}
#'
#' @examples
#' data(Ex4belly)
#' library(survival)
#' km1 <- with(Ex4belly,survfit(Surv(month,evntd)~trt))
#' km1
#' plot(km1)
#' with(subset(Ex4belly,trt==1),pwexpfit(Srv=Surv(month,evntd),array(6,3)))
#' with(subset(Ex4belly,trt==0),pwexpfit(Srv=Surv(month,evntd),array(6,3)))
"Ex4belly"