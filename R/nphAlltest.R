#' nphAlltest A Wrapper Function
#'
#' A wrapper function to run all candidate tests on given data
#'
#' The function runs the following tests then consolidates key statistics into a summary table. \cr
#' Log-Rank test/cox proportional hazard model, max-combo test by Larry, Max-WLR Test by Satrajit and restricted mean test.
#'
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment variable

#' @return The function returns a list with the following components
#' \describe{
#'  \item{lr}{raw output from log-Rank test} 
#'  \item{cox}{raw output from cox PH model}
#'  \item{comboLR}{raw output from max-combo test by Larry}
#'  \item{maxWLR}{raw output from max-WLR test by Satrajit}
#'  \item{rmst}{raw output from restricted mean test}
#'  \item{summary.stat}{concatenated summary statistics including 1-sided p-value, point estimate (mean difference/ratio for RMST tests and hazard ratios for the rest tests), lower and upper bound of 95 percent confidence interval}
#'  } 
#' @examples
#' library(survival)
#' library(survMisc)
#' set.seed(1)
#' indf <- data.frame(trt=rep(c('control','expermental'),100),aval=rexp(200,c(0.002,0.0018)),
#'                    event=rbinom(200,1,0.9))
#' out <- nphAlltest(survival=indf$aval,cnsr=1-indf$event,trt=indf$trt=='expermental')
#' out
#'
#' @export
#' @import data.table survival 
#' 

nphAlltest <- function(survival, cnsr, trt){
  # log-rank test and cox ph model
  lr <- survdiff(Surv(survival,1-cnsr) ~ trt)
  lrstat <- ifelse(diff(lr$obs-lr$exp)>0,-sqrt(lr$chisq),sqrt(lr$chisq))
  cox.fit  <- coxph(Surv(survival,1-cnsr) ~ trt)
  stat1 <- data.frame(Test='Log-Rank',t(c(1-pnorm(lrstat),summary(cox.fit)$conf.int[c(1,3,4)])))
  
  # weighted log rank test (FH)
  rgs <- list(c(0, 0), c(0, 1), c(0, 2), c(0, 3), c(1, 0), c(2, 0), c(2, 2))
  wlr <- fit.combo.wlr(rgs = rgs, time = survival, delta = 1-cnsr, z = trt, draws = 1000,
                       plot.rg = FALSE, print.results = TRUE, outfile = NULL)
  
  stat2 <- data.frame(Test=paste0('Weighted Log-Rank (FH) wt=(',wlr$WLR[,'rho'],',',wlr$WLR[,'gamma'],')'),wlr$WLR[,4:7])
  
  # combination test
  CoxMax <- merge(data.frame(wlr$CoxMax)[,c(1,2,4,5)],wlr$WLR[,c('rho','gamma','hrL.band','hrU.band')],by=c('rho','gamma'))
  stat3 <- data.frame(Test=paste0('Max-Combination Test (FH) wt=(',CoxMax[,'rho'],',',CoxMax[,'gamma'],')'),CoxMax[1,3:6])
  
  # max-WLR test
  comboLR <- rm.combo.WLRmax(time=survival,status=1-cnsr,arm=trt,wt=list(a1=rgs[[2]],a2=rgs[[3]],a3=rgs[[4]],a4=rgs[[5]],a5=rgs[[6]],a6=rgs[[7]]),max=T)
  combstat <- unlist(comboLR)
  stat4 <- data.frame(Test=paste0('Max-WLR Test (FH) wt=(',rgs[[combstat[5]]][1],',',rgs[[combstat[5]]][2],')'),t(combstat[c(1,2,6,7)]))
  
  # rmst and rmtl
  rmst <- rmst2(survival,1-cnsr,trt)
  stat5<- data.frame(Test=c('RMST Difference','RMST Ratio','RMTL Ratio'),rmst$unadjusted.result[,c(4,1:3)],row.names=NULL)
  if (stat5[1,3]>0) stat5[,2] <-stat5[,2]/2 else stat5[,2] <- 1-stat5[,2]/2  # convert to 1-sided
  col <- c('Test','p-Value','Estimate','Lower','Upper')
  
  out.stat <- rbind(setNames(stat1,col),setNames(stat2,col),setNames(stat3,col),setNames(stat4,col),setNames(stat5,col))
  return (list(lr=lr,cox=cox.fit,maxWLR=wlr,comboLR=comboLR,rmst=rmst,summary.stat=out.stat))
}





