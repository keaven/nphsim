## ----setup, include=FALSE------------------------------------------------
options(scipen=999)
knitr::opts_chunk$set(echo = TRUE)
library(nphsim)
library(survminer)
library(survival)
library(knitr)

## ----fig.height=5,fig.width=6--------------------------------------------
fit1 <- survfit(Surv(month, evntd) ~ trt,
               data = Ex1delayedEffect)
fit2 <- survfit(Surv(month, evntd) ~ trt,
               data = Ex2delayedEffect)
Ex1KMPlot<-ggsurvplot(fit1, data = Ex1delayedEffect, risk.table = TRUE, break.time.by=3)
Ex2KMPlot<-ggsurvplot(fit2, data = Ex2delayedEffect, risk.table = TRUE, break.time.by=3)
Ex1KMPlot
Ex2KMPlot

## ------------------------------------------------------------------------
kable(Ex1Rate1 <- with(subset(Ex1delayedEffect,trt==1), pwexpfit(Surv(month,evntd))),
      caption="Test case 1 piecewise exponential fit, experimental group")

## ----echo=FALSE----------------------------------------------------------
kable(Ex1Rate0 <- with(subset(Ex1delayedEffect,trt==0), pwexpfit(Surv(month,evntd))),
      caption="Test case 1 piecewise exponential fit, control group")

## ------------------------------------------------------------------------
Ex1RateAll <-  with(Ex1delayedEffect, pwexpfit(Surv(month,evntd)))
pchisq(sum(Ex1RateAll$m2ll)-sum(Ex1Rate0$m2ll+Ex1Rate1$m2ll),df=4,lower.tail=FALSE)
lr1 <- with(Ex1delayedEffect, survdiff(formula = Surv(month, evntd) ~ trt))$chisq
pchisq(lr1,df=1,lower.tail=FALSE)

## ----echo=FALSE----------------------------------------------------------
Ex2Rate1 <- with(subset(Ex2delayedEffect,trt==1), pwexpfit(Surv(month,evntd)))
Ex2Rate0 <- with(subset(Ex2delayedEffect,trt==0), pwexpfit(Surv(month,evntd)))
rate2all <-  with(Ex2delayedEffect, pwexpfit(Surv(month,evntd)))
p2df4 <- pchisq(sum(rate2all$m2ll)-sum(Ex2Rate0$m2ll+Ex2Rate1$m2ll),df=4,lower.tail=FALSE)
lr2 <- with(Ex1delayedEffect, survdiff(formula = Surv(month, evntd) ~ trt))$chisq
p2lr <- pchisq(lr2, df=1, lower.tail=FALSE)

## ----echo=FALSE----------------------------------------------------------
kable(Ex2Rate1, caption="Test case 2 piecewise exponential fit, experimental group")
kable(Ex2Rate0, caption="Test case 2 piecewise exponential fit, control group")

## ------------------------------------------------------------------------
# example 1 piecewise hazard ratios
Ex1Rate1$rate / Ex1Rate0$rate
# example 2 piecewise hazard ratios
Ex2Rate1$rate / Ex2Rate0$rate


## ----echo=FALSE----------------------------------------------------------
# common rate for first 3 months
Ex1RateAllMonth3 <- Ex1RateAll[1,]
Ex1RateAllMonth3$Treatment <- 'Combined'
Ex1RateAllMonth3$Period <- '< 3 months'
# Combine time periods after 3 months into 1 for experimental treatment
Ex1Rate1PostMonth3 <- with(subset(Ex1delayedEffect,trt==1),
                           pwexpfit(Surv(month,evntd),intervals=c(3,Inf)))[2,]
Ex1Rate1PostMonth3$Treatment <- 'Experimental'
Ex1Rate1PostMonth3$Period <- '> month 3'
# Repeat for control
Ex1Rate0PostMonth3 <- with(subset(Ex1delayedEffect,trt==0),
                           pwexpfit(Surv(month,evntd),intervals=c(3,Inf)))[2,]
Ex1Rate0PostMonth3$Treatment <- 'Control'
Ex1Rate0PostMonth3$Period <- '> month 3'
# combine and print table
Ex1RateSimple <- rbind(Ex1RateAllMonth3,Ex1Rate0PostMonth3,Ex1Rate1PostMonth3)
kable(Ex1RateSimple[,-1], caption="2-piece exponential fit with HR=1 for first 3 months, test case 1.")

## ----echo=FALSE----------------------------------------------------------
Ex1PWEMonth3 <- exp(-Ex1RateSimple$rate[1]*(0:30)*.1)
Ex1PWE1PostMonth3 <- Ex1PWEMonth3[31]*exp(-Ex1RateSimple$rate[3]*(1:120)*.1)
Ex1PWE0PostMonth3 <- Ex1PWEMonth3[31]*exp(-Ex1RateSimple$rate[2]*(1:120)*.1)
Ex1PWESurvData <- rbind(data.frame(time=(0:150)*.1,surv=c(Ex1PWEMonth3,Ex1PWE0PostMonth3),strata='trt=0'),
                        data.frame(time=(0:150)*.1,surv=c(Ex1PWEMonth3,Ex1PWE1PostMonth3),strata='trt=1'))
Ex1KMPlot$plot+geom_line(data=Ex1PWESurvData,aes(x=time,y=surv,col=strata))

