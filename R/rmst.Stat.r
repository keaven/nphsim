#' Restricted Mean Survival Test Statistic
#'
#' RMST method with user specified trucation time.
#'
#' Restricted mean survival test at specified cutoff time. This is adopted from \code{\link[survRM2:rmst2]{rmst2}}. 
#'
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment varaible
#' @param stra stratification variable. Default is \code{NULL}
#' @param fparam the cutoff time for RMST analysis. If larger than the minimum of the largest observed time 
#' on each of the two arms, the minimum value will be used.
#' @return The function return a list with the follow components
#' \describe{
#'  \item{pval}{One-sided p-Value from RMST comparing treatment arm to control arm} 
#'  \item{tau}{time cutoff used for each analysis}
#'  \item{est}{Estimated difference between treatment and control arm}
#'  \item{estlb, estub}{Lower and upper bound of the estimated difference}
#'  }
#'   
#' @examples
#' # RMST on the simulated data
#' library(survRM2)
#' medC = 6 
#' hr <- 0.7
#' intervals <- NULL 
#' gamma <- c(2.5, 5,  7.5,  10) ## a ramp-up enrollment
#' R     <- c(2  , 2,  2  ,  6 ) ## enrollment period: total of 12 months
#' eta <- -log(0.99) ## 1% monthly dropout rate
#' sim1 <- nphsim(nsim=1,lambdaC=log(2)/medC,lambdaE=log(2)/medC*hr, ssC=300,ssE=300,intervals=intervals,gamma=gamma, R=R,eta=eta)
#' test1 <- simtest(x=sim1, anaD=c(250,300), method=rmst.Stat,fparam=12)
#' test1$result[]
#' 
#' # direct function call (without cutoff)
#' rmst.Stat(surv=sim1$simd$survival, cnsr=sim1$simd$cnsr, trt=sim1$simd$treatment,fparam=12)
#' 
#' @import survRM2
#' @export
rmst.Stat <- function(survival, cnsr, trt, stra = NULL, fparam = NULL) {
  idx = arm == 0
  tt = time[idx]
  tau0max = max(tt)
  idx = arm == 1
  tt = time[idx]
  tau1max = max(tt)
  tau_max = min(tau0max, tau1max)
  tau <- ifelse(tau_max < fparam, tau_max, fparam)
  a <- rmst2(time = survival, status = 1 - cnsr, arm = (trt == "experiment"), tau = tau)
  b <- a$unadjusted.result
  pval <- ifelse(sign(b[1, 1]), b[1, 4]/2, (1 - b[1, 4])/2)
  y <- list(pval = pval, tau = tau, est = b[1, 1], estlb = b[1, 2], estub = b[1, 3])
  return(y)
}

