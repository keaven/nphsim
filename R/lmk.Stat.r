#' Landmark analysis
#'
#' Landmark (survival curve comparison at fixed point in time) analysis based on Klein 2007 (statistics in Medicine)
#'
#' The p-Values from the following transformations of survival are included:
#' \describe{
#'  \item{naive}{No transformation} 
#'  \item{log}{\code{log(S)}}
#'  \item{loglog}{\code{log(-log(S))}}
#'  \item{arcsin}{\code{arcsin(sqrt(S))}}
#'  \item{logit}{\code{log(S/(1-S))}}
#'  }
#'
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment varaible. Accepted values are either "experiment" or "control"
#' @param stra stratification variable. Default is \code{NULL} (currently not implemented)
#' @param fparam a list input
#' \describe{
#'  \item{\code{fparam$lmk}}{time at which the landmark analysis is done} 
#'  \item{\code{fparam$lmktype}}{test name of which the p-Value will be stored in the variable \code{pval}. 
#'  Accepted values are "naive", "log", "loglog", "arcsin", "logit"}
#'  } 
#' @return The function return a list with the follow components
#' \describe{
#'  \item{pval}{One-sided p-Value from user specified test} 
#'  \item{pval_[XXX]}{One-sided test from various transformations}
#'  }
#'   
#' @examples
#' # Lankmark analysis on the simulated data
#' library(survMisc)
#' medC = 6 
#' hr <- c(1, 0.6)
#' intervals <- 3 
#' gamma <- c(2.5, 5,  7.5,  10) ## a ramp-up enrollment
#' R     <- c(2  , 2,  2  ,  6 ) ## enrollment period: total of 12 months
#' eta <- -log(0.99) ## 1% monthly dropout rate
#' sim1 <- nphsim(nsim=1,lambdaC=log(2)/medC,lambdaE=log(2)/medC*hr, ssC=300,ssE=300, 
#'                intervals=intervals,gamma=gamma, R=R,eta=eta)
#' test1 <- simtest(x=sim1, anaD=c(250,300), method=lmk.Stat,fparam=list(lmk=9, lmktype='loglog'))
#' test1$result[]
#' 
#' # direct function call (without cutoff)
#' lmk.Stat(surv=sim1$simd$survival, cnsr=sim1$simd$cnsr, trt=sim1$simd$treatment, 
#'          fparam=list(lmk=9, lmktype='loglog'))
#' 
#' @export
#' @import data.table
#' 
lmk.Stat <- function (survival, cnsr, trt, stra = NULL, fparam) {
  d <- data.table(survival = survival, cnsr = cnsr, trt = trt)
  x <- ten(survfit(Surv(survival, 1 - cnsr) ~ trt, data = d))
  s <- sf(x)
  ts <- data.table(s[t < fparam$lmk][, .SD[.N], by = cg])
  ts[, sigma2 := ifelse(S != 0, Sv / (S ^ 2), 1e5)]
  
  n1 <- c("naive", "log", "loglog", "arcsin", "logit")
  n2 <- c("transform", "pChisq", "chiSq")
  
  res1 <- data.table(ts, transform = rep(n1, each = 2))
  
  res1[transform == "loglog", a1 := log(-log(S))][, a2 := sigma2 / log(S) ^ 2]
  res1[transform == "naive", a1 := S][, a2 := sigma2 * S ^ 2]
  res1[transform == "log", a1 := log(S)][, a2 := sigma2]
  res1[transform == "arcsin", a1 := asin(sqrt(S))][, a2 := sigma2 * S /
                                                     (4 * (1 - S))]
  res1[transform == "logit", a1 := log(S / (1 - S))][, a2 := sigma2 / (1 - S) ^ 2]
  res <- res1[, .(chisq = sum(a1 * c(1, -1)) ^ 2 / sum(a2 * c(1, 1))), by = transform]
  res[, pval := (1 - pchisq(chisq, 1)) / 2]
  pvals <- transpose(res[, .(pval)])
  setnames(pvals, paste0("pval_", n1))
  pval = res[transform == fparam$lmktype, pval]
  y <- cbind(pval, pvals)
  
  return(y)
  
}

