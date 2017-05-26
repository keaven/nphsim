#' LR Test Statistic
#'
#' Tests if there is a difference between two survival curves using logrank test. One sided p-value, Z statistics, HR from cox model and standard error
#' of the HR are generated.
#'
#' test.lr is used when method='LR'
#'
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment varaible. Accepted values are either "experiment" or "control"
#' @param stra stratification variable. Default is \code{NULL}
#' @return The function return a list with the follow components
#' \describe{
#'  \item{pval}{One-side p-Value from the logrank test} 
#'  \item{z}{z statistics from the logrank test}
#'  \item{hr}{Hazard ratio from the cox proportional hazard model}
#'  \item{sehr}{Standard error of the hazard ratio}
#'  }
#'   
#' @examples
#' # logrank test on the simulated data
#' library(survival)
#' medC = 6 
#' hr <- 0.7
#' intervals <- NULL 
#' gamma <- c(2.5, 5,  7.5,  10) ## a ramp-up enrollment
#' R     <- c(2  , 2,  2  ,  6 ) ## enrollment period: total of 12 months
#' eta <- -log(0.99) ## 1% monthly dropout rate
#' sim1 <- nphsim(nsim=1,lambdaC=log(2)/medC,lambdaE=log(2)/medC*hr, ssC=300,ssE=300,intervals=intervals,gamma=gamma, R=R,eta=eta)
#' test1 <- simtest(x=sim1, anaD=300, method='LR')
#' test1$result[]
#' 
#' # direct function call (without cutoff)
#' test.lr(surv=sim1$simd$survival, cnsr=sim1$simd$cnsr, trt=sim1$simd$treatment)
#' 
#' 
#' @import survival
#' @export
test.lr <- function(surv, cnsr, trt, stra = NULL) {
  ## logrank test p-value
  if (is.null(stra)) {
    lr <- survdiff(Surv(surv, 1 - cnsr) ~ trt)
    cox <- coxph(Surv(surv, 1 - cnsr) ~ trt)
  } else {
    lr <- survdiff(Surv(surv, 1 - cnsr) ~ trt + strata(stra))
    cox <- coxph(Surv(surv, 1 - cnsr) ~ trt + strata(stra))
  }
  trt_better <- sign(lr$obs[2] - lr$exp[2])
  pval <- ifelse(sign(trt_better), (1 - pchisq(lr$chisq, length(unique(trt)) - 1))/2, pchisq(lr$chisq, length(unique(trt)) - 1)/2)
  y <- list(pval = pval, z = sign(lr$obs[2] - lr$exp[2]) * sqrt(lr$chisq), hr = exp(cox$coefficients), sehr = sqrt(diag(cox$var)))
  return(y)
}
