#' nphAlltest A Wrapper Function
#'
#' A wrapper function to run all candidate tests on given data
#'
#' The function runs the following tests then consolidates key statistics into a summary table. \cr
#' Log-Rank test/cox proportional hazard model, max-combo test by Larry, Max-WLR Test by 
#' Satrajit and restricted mean test.
#'
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment variable

#' @return The function returns a list with the following components
#' \describe{
#'  \item{lr}{raw output from log-Rank test} 
#'  \item{cox}{raw output from cox PH model}
#'  \item{maxWLR}{raw output from max-combo test by Larry}
#'  \item{comboLR}{raw output from max-WLR test by Satrajit}
#'  \item{rmst}{raw output from restricted mean test}
#'  \item{summary.stat}{concatenated summary statistics including 1-sided p-value, point estimate 
#'  (mean difference/ratio for RMST tests and hazard ratios for the rest tests), lower and upper 
#'  bound of 95 percent confidence interval}
#'  } 
#' @examples
#' library(survival)
#' library(survMisc)
#' library(mvtnorm)
#' library(Matrix)
#' library(dplyr)
#' set.seed(1)
#' indf <- data.frame(trt=rep(c('control','experimental'),100),aval=rexp(200,c(0.002,0.0018)),
#'                    event=rbinom(200,1,0.9))
#' out <- nphAlltest(survival=indf$aval,cnsr=1-indf$event,trt=indf$trt)
#' out$summary.stat
#' out <- nphAlltest(survival=indf$aval,cnsr=1-indf$event,trt=indf$trt,rgs=list(c(0,0),c(0,1)))
#' out
#' with(Ex1delayedEffect, nphAlltest(month,1-evntd,trt=factor(trt,labels=c("control","experimental"))))
#'
#' @export
#' @import data.table survival 
#' 

nphAlltest <- function(survival,cnsr,trt,rgs=list(c(0,0),c(0,1),c(0,2),c(0,3),c(1,0),c(2,0),c(2,2))){
  # log-rank test and cox ph model
  lr <- unlist(test.lr(survival, cnsr, trt))
  lr.ci <- exp(c(log(lr[3])-qnorm(0.975)*lr[4],log(lr[3])+qnorm(0.975)*lr[4]))
  stat1 <- data.frame(Test='Log-Rank',t(c(lr[c(1,3)],lr.ci)))
  
  # weighted log rank test (FH)
  wlr <- fit.combo.wlr(rgs = rgs, time = survival, delta = 1-cnsr, z = trt=="experimental", draws = 1000,
                       plot.rg = FALSE, print.results = TRUE, outfile = NULL)
  stat2 <- data.frame(Test=paste0('Weighted Log-Rank (FH) wt=(',wlr$WLR[,'rho'],',',
                                  wlr$WLR[,'gamma'],')'),wlr$WLR[,4:7])
  
  # max combination test
  CoxMax <- merge(data.frame(wlr$CoxMax)[,c(1,2,4,5)],wlr$WLR[,c('rho','gamma','hrL.band','hrU.band')],
                  by=c('rho','gamma'))
  stat3 <- data.frame(Test=paste0('Max-Combination Test (FH) wt=(',CoxMax[,'rho'],',',
                                  CoxMax[,'gamma'],')'),CoxMax[1,3:6])
  
  # Asymptotic combination test
  AsympCombo <- rm.combo.WLRmax(time=survival,status=1-cnsr,arm=(trt=='experimental'), wt=rgs[-1], 
                                adjust.methods="asymp",one.sided=T)
  combstat <- unlist(AsympCombo[c('pval','hr.est','hr.low.adjusted','hr.up.adjusted','max.index')])
  stat4 <- data.frame(Test=paste0('Asymp-Combination Test (FH) wt=(',rgs[[combstat[5]]][1],',',
                                  rgs[[combstat[5]]][2],')'),t(combstat[1:4]))
  
  
  # restricted mean difference
  rmst <- unlist(rmst.Stat(survival,cnsr,trt,fparam=max(survival[cnsr==0])))
  stat5 <- data.frame(Test='RMST Difference',t(rmst[c(1,3:5)]),row.names=NULL)
  
  col <- c('Test','p-Value','Estimate','Lower','Upper')
  out.stat <- rbind(setNames(stat1,col),setNames(stat2,col),setNames(stat3,col),
                    setNames(stat4,col),setNames(stat5,col))
  return (list(MaxCombo=wlr,AsympCombo=AsympCombo,summary.stat=out.stat))
}




