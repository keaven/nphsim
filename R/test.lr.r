#' LR Test Statistic
#'
#' Description here
#'
#' test.lr is used when method='LR'
#'
#' @param survival parameter description
#' @param cnsr parameter description
#' @param trt  parameter description
#' @param stra parameter description
#' @examples
#' # TBD
#' @import survival
#' @export
test.lr<-function (surv,cnsr,trt,stra=NULL) {
  ## logrank test p-value
  if (is.null(stra)){
    lr <- survdiff(Surv(surv,1-cnsr)~trt)
    cox <- coxph(Surv(surv,1-cnsr)~trt)
  } else{
    lr <- survdiff(Surv(surv,1-cnsr)~trt + strata(stra))
    cox <- coxph(Surv(surv,1-cnsr)~trt+strata(stra))
  }
  pval <- (1 - pchisq(lr$chisq, length(unique(trt))-1))/2
  y<-list(pval=pval, hr=exp(cox$coefficients), sehr=sqrt(diag(cox$var)))
  return(y)
}
