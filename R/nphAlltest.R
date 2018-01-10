#' nphAlltest A Wrapper Function
#'
#' A wrapper function to run all candidate tests for a treatment difference on given data
#'
#' The function runs the following tests then consolidates key statistics into a summary table. \cr
#' Log-Rank test/cox proportional hazard model, maximum weighted logrank (asymptotic p-value computation), 
#' restricted mean survival difference and ratio as well as restricted mean time lost ratio; note that the rmst2 
#' default cutoff time of the minimum of the maximum observed time in each group.
#'
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment variable
#' @param cutpt a cutpoint or cutpoints to use in defining time intervals for piecewise hazard ratio estimation
#' @param KMtitle title for Kaplan-Meier plot
#' @param cumhaztitle title for cumulative hazard plot
#' @return The function returns a list with the following components
#' \describe{
#'  \item{summary.stat}{concatenated summary statistics including 1-sided p-value, point estimate 
#'  (mean difference/ratio for RMST tests and hazard ratios for the rest tests), lower and upper 
#'  bound of 95 percent confidence interval}
#'  \item{hr.pe}{Piecewise Cox model estimate, SE and 95 percent confidence limits for
#'  hazard ratio by intervals specified in cutpt}
#'  \item{km.plot}{Kaplan-Meier plot}
#'  \item{cumhaz.plot}{Cumulative hazard plot}
#'  } 
#' @examples
#' library(survival)
#' library(survMisc)
#' library(mvtnorm)
#' library(Matrix)
#' library(dplyr)
#' library(survminer)
#' library(survRM2)
#' set.seed(1)
#' indf <- data.frame(trt=rep(c('control','experimental'),100),aval=rexp(200,c(0.002,0.0018)),
#'                    event=rbinom(200,1,0.9))
#' out <- nphAlltest(survival=indf$aval,cnsr=1-indf$event,trt=indf$trt)
#' out$summary.stat
#' out$km.plot
#' out$cumhaz.plot
#' x1out <- with(Ex1delayedEffect, nphAlltest(month,1-evntd,trt=factor(trt,labels=c("control","experimental"))))
#' x1out$summary.stat
#' @export
#' @import data.table survival 
#' 

nphAlltest <- function(survival,cnsr,trt,cutpt=median(survival),
                       KMtitle="Kaplan-Meier Plot",
                       cumhaztitle="Cumulative Hazard Plot"){
  # log-rank test and cox ph model
  lr <- unlist(test.lr(survival, cnsr, trt))
  lr.ci <- exp(c(log(lr[3])-qnorm(0.975)*lr[4],log(lr[3])+qnorm(0.975)*lr[4]))
  stat1 <- data.frame(Test='Log-Rank',t(c(lr[c(1,3)],lr.ci)))
  
  # weighted log rank test (FH)
  wlr <- function (wt) c(wt[1],wt[2],unlist(rm.combo.WLRmax(time=survival,status=1-cnsr,
                                                            arm=(trt=='experimental'),
                                                            wt=list(a1=c(wt[1],wt[2])), max=F)))
  wlr.est <- data.frame(t(apply(matrix(c(0,1,0,2,0,3,1,0,2,0,1,1,2,2),nrow=2),2,wlr)))
  stat2 <- data.frame(Test=paste0('Weighted Log-Rank (FH) wt=(',wlr.est[,1],',',wlr.est[,2],')'),wlr.est[,3:6])
  
  # Asymptotic combination test with (0,0),(0,1)
  rgs1 <- list(c(0,0),c(0,1))
  AsympCombo1 <- rm.combo.WLRmax(time=survival,status=1-cnsr,arm=(trt=='experimental'), wt=rgs1[-1], 
                                 adjust.methods="asymp",one.sided=T)
  combstat1 <- unlist(AsympCombo1[c('pval','hr.est','hr.low.adjusted','hr.up.adjusted','max.index')])
  stat3 <- data.frame(Test=paste0('FH Combination ((0,0),(0,1)) wt=(',rgs1[[combstat1[5]]][1],',',
                                  rgs1[[combstat1[5]]][2],')'),t(combstat1[1:4]))
  
  # Asymptotic combination test with (0,0),(0,1),(1,0),(1,1)
  rgs2 <- list(c(0,0),c(0,1),c(1,0),c(1,1))
  AsympCombo2 <- rm.combo.WLRmax(time=survival,status=1-cnsr,arm=(trt=='experimental'), wt=rgs2[-1], 
                                 adjust.methods="asymp",one.sided=T)
  combstat2 <- unlist(AsympCombo2[c('pval','hr.est','hr.low.adjusted','hr.up.adjusted','max.index')])
  stat4 <- data.frame(Test=paste0('FH Combination ((0,0),(0,1),(1,0),(1,1)) wt=(',rgs2[[combstat2[5]]][1],
                                  ',',rgs2[[combstat2[5]]][2],')'),t(combstat2[1:4])) 
  
  # restricted mean difference
  rmst <- rmst2(survival,1-cnsr,trt=='experimental')
  stat5<- data.frame(Test=c('RMST Difference','RMST Ratio','RMTL Ratio'),rmst$unadjusted.result[,c(4,1:3)],row.names=NULL)
  if (stat5[1,3]>0) stat5[,2] <-stat5[,2]/2 else stat5[,2] <- 1-stat5[,2]/2  # convert to 1-sided
  
  
  col <- c('Test','p-Value','Estimate','Lower','Upper')
  out.stat <- rbind(setNames(stat1,col),setNames(stat2,col),setNames(stat3,col),setNames(stat4,col),setNames(stat5,col))
  
  # piecewise HR
  hr.pe <- cox_pw(time=survival,status=1-cnsr,arm=ifelse(trt=='experimental',1,0), cutpt)$est
  
#   plot cumulative hazard and survival curves
  d <- data.frame(survival=survival,status=1-cnsr,trt=trt)
  fit <- survfit(Surv(survival, status)~trt,data=d)
  cumhaz.plot <- ggsurvplot(fit,data=d,
                            title=cumhaztitle,
                            fun="cumhaz",
                            palette=c("#1B9E77","#D95F02"),
                            risk.table=TRUE)
  km.plot <- autoplot(fit,
                      title=KMtitle,
                      xLab='Time',
                      yLab='Survival Probability',
                      tabTitleSize=10, 
                      nRiskSize=4,
                      yScale='frac',
                      censSize=2)
  
  return (list(summary.stat=out.stat,hr.pe=hr.pe,cumhaz.plot=cumhaz.plot,
               km.plot=km.plot))
}
