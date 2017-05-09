sapply(c(
  #'doParallel',
  #'plyr',
  #'reshape2',
  'survival',
  'survRM2',
  'survMisc'
), require, character.only=TRUE)

#######################################################################################
## test.lr is used when method='LR'
#######################################################################################
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

#######################################################################################
## Example of User provided testing function
## RMST method with user specified trucation time
#######################################################################################
rmst.Stat<-function (survival,cnsr,trt,stra=NULL,fparam=NULL) {
  a<-rmst2(time=survival,status=1-cnsr,arm=(trt=='experiment'),tau=fparam$tt)
  b<-a$unadjusted.result
  pval <- b[1,4]/2
  y<-list(pval=pval, est=b[1,1], estlb=b[1,2],estub=b[1,3])
  return(y)
}


###################################################################################
##                                                                               ##
## Produce the Weighted Kaplan-Meier test statistic(borrowed from ELSTIC NPH)    ##                                                                               ##
## Input : (1) dataframe 'indata' contains the following variables:              ##
##             Event.T -> Event time (failure/censor)                            ##
##              Status -> 0 = Censored, 1 = Event of Interest                    ##
##                   Z -> 0 = Control, 1 = Active                                ##
##                                                                               ##
## Dependence:{survival}                                                         ##
##                                                                               ##
## Note: (1) Test statistic is constructed based on the paper of Pepe and        ##
##           Fleming (1989).                                                     ##         
##       (2) Code returns the Z statistic and p-value.                           ##
##       (3) Code is adapted from '6test.core.r'.                                ##
##       (4) Although it is of interest to if the statistic will reduce to RMST  ##
##           when the weight function = 1, the stability conditions proposed by  ##
##           Pepe and Fleming (1989) rule out a constant weight function when    ##
##           censoring is present.                                               ##
##                                                                               ##
###################################################################################
WKM.Stat <- function(survival,cnsr,trt,stra=NULL,fparam=NULL){
  indata<-data.table(Event.T=survival,Status=1-cnsr,Z=(trt=='experiment'))
  data.g1 <- indata[Z==1, ]
  data.g2 <- indata[Z==0, ]
  n1 <- data.g1[,.N]
  n2 <- data.g2[,.N]
  n <- n1 + n2
  
  
  ## Generate a survfit object based on the failure times for each group
  
  sur.fit.g1 <- summary(survfit(Surv(Event.T, Status == 1) ~ 1, data = data.g1))
  sur.fit.g2 <- summary(survfit(Surv(Event.T, Status == 1) ~ 1, data = data.g2))
  
  ## Compute the Kaplan-Meier estimate of the survivor function at all event times
  ##   (including censoring time). Use linear interpolate function 'approxfun' to 
  ##   create the KM.sx() function with option 'f=0' so that KM.sx() is a  
  ##   Right-continuous with Left limit (Cadlag) function.
  
  KM.s1.f <- approxfun(sur.fit.g1$time, sur.fit.g1$surv, method = "constant", 
                       yleft = 1, rule = 2, f = 0)
  KM.s2.f <- approxfun(sur.fit.g2$time, sur.fit.g2$surv, method = "constant", 
                       yleft = 1, rule = 2, f = 0)
  
  
  ## Generate a grid for calculating the integrals. Since the integrand is a Step 
  ##   function with its grid (or discontinuous) points being  continuous from the 
  ##   Left or Right, it is better to use the midpoint sum to approximate the integral.
  
  event.time <- sort(indata$Event.T)
  
  integrand.grid <- as.vector(
    apply( matrix(cbind(c(0, event.time[-n]), event.time), 
                  ncol = 2, byrow = F),
           1, mean
    )
  )
  
  
  ## K-M estimates' values at midpoints between event times (including the failure 
  ##   times and censoring times)
  
  KM.s1 <- KM.s1.f(integrand.grid)
  KM.s2 <- KM.s2.f(integrand.grid)
  
  
  
  ## ----------------------- calculate the weight function ------------------------- ##
  ##                          Section 4 in Pepe & Fleming                            ##
  
  cen.fit.g1 <- summary(survfit(Surv(Event.T, Status == 0) ~ 1, data = data.g1))
  cen.fit.g2 <- summary(survfit(Surv(Event.T, Status == 0) ~ 1, data = data.g2))
  
  ## Compute the Kaplan-Meier estimate of the censoring probability at all event times.
  ##   (including (s1, s2)). Use linear interpolate function 'approxfun' to create the 
  ##   KM.cx() function. Since the survival functions for censoring are only used in 
  ##   the Weight function where the weight function uses the Left-continuous version 
  ##   (cf. Section 4), the code will use the option 'f=1' to make KM.cx() to be 
  ##   Left-continuous  with Right limit (Caglad) function.
  
  n1.cen <- length(cen.fit.g1$time)
  n2.cen <- length(cen.fit.g2$time)
  
  if ( n1.cen >= 1 ) {
    KM.cm1.f <- approxfun( c(cen.fit.g1$time, max(cen.fit.g1$time) + 1), 
                           c(1, cen.fit.g1$surv), method = "constant", 
                           yleft = 1, yright = min(cen.fit.g1$surv), rule = 2, f = 1)
  } else {
    KM.cm1.f <- approxfun( c(0, max(event.time)), c(1, 1), method = "constant", 
                           yleft = 1, rule = 2, f = 1)
  }
  
  if ( n2.cen >= 1 ) {
    KM.cm2.f <- approxfun( c(cen.fit.g2$time, max(cen.fit.g2$time) + 1), 
                           c(1, cen.fit.g2$surv), method = "constant", 
                           yleft = 1, yright = min(cen.fit.g2$surv), rule = 2, f = 1)
  } else {
    KM.cm2.f <- approxfun( c(0, max(event.time)), c(1, 1), method = "constant", 
                           yleft = 1, rule = 2, f = 1)
  }  
  
  KM.cm1 <- KM.cm1.f(integrand.grid)
  KM.cm2 <- KM.cm2.f(integrand.grid)
  
  weight <- ifelse(KM.cm1 + KM.cm2 == 0, 0, (n * KM.cm1 * KM.cm2)/(n1 * KM.cm1 + n2 * KM.cm2))
  
  
  
  ## ----------------------- calculate the Test Statistic ------------------------- ##
  ##                        Formula (3.3) in Pepe & Fleming                         ##
  
  ############################################################
  ## ---------  Numerator of the Score statistic  --------- ##
  ##            Equation (3.3) in Pepe & Fleming            ##
  ############################################################
  
  num <- sqrt((n1 * n2)/n) * 
    sum(weight * (KM.s1 - KM.s2) * diff(c(0, event.time)), na.rm = T)
  
  
  ## Denominator(standard deviation) of the Score statistic using the 
  ##     pooled variance estimate
  
  KM.pool <- summary(survfit(Surv(Event.T, Status) ~ 1, data = indata))
  
  
  KM.S.f <- approxfun( KM.pool$time, KM.pool$surv, method = "constant", 
                       yleft = 1, rule = 2, f = 0)
  KM.Sm.f <- approxfun( c(KM.pool$time, max(KM.pool$time) + 1), c(1, KM.pool$surv), 
                        method = "constant", yleft = 1, yright = min(KM.pool$surv),
                        rule = 2, f = 1)
  
  KM.S <- KM.S.f(integrand.grid)
  KM.Sm <- KM.Sm.f(integrand.grid)
  
  A.seq <- cumsum(diff(c(0, event.time)) * weight * KM.S)
  A <- (A.seq[n] - A.seq)
  
  
  #       delta <- indata$Status[order(indata$Event.T)] == 1
  #
  #       variance <- - sum( A[delta]^2 / (c(1, KM.pool$surv[-length(KM.pool$surv)]))^2 
  #                          * (weight[delta])^(-1) 
  #                          * diff(c(1, KM.pool$surv)), 
  #                         na.rm = T
  #                     )
  
  
  variance <- - sum( A^2 / (KM.S * KM.Sm) 
                     * (weight)^(-1)
                     * diff(c(KM.Sm, min(KM.pool$surv))), 
                     na.rm = T
  )
  
  denominator <- sqrt(variance)
  
  stat <- num / denominator
  
  return(list(
    Stat = stat,
    Df = NA,
    pval = (1 - pnorm(stat))
  ))
}


#######################################################################################
## Weighted logrank adopted from survMisc package
#######################################################################################
wlr.Stat<-function (survival,cnsr,trt,stra=NULL,fparam) {
  d<-data.table(survival=survival,cnsr=cnsr,trt=trt)
  x<-ten(survfit(Surv(survival, 1-cnsr) ~ trt,data=d))
  t1 <- x[e>0, t, by=t][, t]
  wt1 <- data.table::data.table(array(data=1, dim=c(length(t1),10)))
  
  n1 <- c("1", "n", "sqrtN", "S1", "S2", "FH01","FH10","FH11","FH","APPLE")
  data.table::setnames(wt1, n1)
  ## Gehan-Breslow generalized Wilcoxon, weight = n
  data.table::set(wt1, j="n",
                  value=x[e>0, max(n), by=t][, V1])
  ## Tarone-Ware, weight = sqrt(n)
  data.table::set(wt1, j="sqrtN",
                  value=wt1[, sqrt(.SD), .SDcols="n"])
  ## Peto-Peto, weight = S(t) = modified estimator of survival function
  data.table::set(wt1, j="S1",
                  value=cumprod(x[e > 0, 1 - sum(e) / (max(n) + 1), by=t][, V1]))
  ## modified Peto-Peto (by Andersen), weight = S(t)n / n+1
  data.table::set(wt1, j="S2",
                  value=wt1[, S1] * x[e > 0, max(n) / (max(n) + 1), by=t][, V1])
  ## Fleming-Harrington
  S3 <- sf(x=x[e>0, sum(e), by=t][, V1], n=x[e>0, max(n), by=t][, V1], what="S")
  S3 <- c(1, S3[seq.int(length(S3) - 1L)])
  wt1[, FH01 := (1 - S3)
      ][, FH10 := S3
      ][, FH11 := S3 * (1 - S3)
      ][, FH :=  S3^fparam$FH[1] * ((1 - S3)^fparam$FH[2])]  ## FH p, q
  
  ## simplified APPLE, w1=0 and w2=1, separate at delayd separation point
  data.table::set(wt1, j="APPLE",
                  value=x[e>0, ifelse(t<fparam$APPLE,0,1), by=t][, V1])
  
  n2 <- c("W", "Q", "Var", "Z",
          "pNorm", "chiSq", "df", "pChisq")
  res1 <- data.table(matrix(0, nrow=ncol(wt1), ncol=length(n2)))
  setnames(res1, n2)
  set(res1, j=1L, value=n1)
  predict(x)
  ## events minus predicted
  eMP1 <- attr(x, "pred")
  eMP1 <- eMP1[rowSums(eMP1) > 0, ]
  ## covariance
  COV(x)
  cov1 <- attr(x, "COV")
  if (is.null(dim(cov1))) {
    cov1 <- cov1[names(cov1) %in% t1]
  } else {
    ## 3rd dimension = times
    cov1 <- cov1[, , dimnames(cov1)[[3]] %in% t1]
  }
  eMP1 <- unlist(eMP1[, .SD, .SDcols=(length(eMP1) - 1L)])
  data.table::set(res1, j="Q",
                  value=colSums(wt1 * eMP1))
  data.table::set(res1, j="Var",
                  value=colSums(wt1^2 * cov1))
  
  res1[, "Z" := Q / sqrt(Var)]
  res1[, "pNorm" := 2 * (1 - stats::pnorm(abs(Z)))]
  res1[, "chiSq" := Q^2 / Var]
  res1[, "df" := 1]
  res1[, "pChisq" := 1 - stats::pchisq(chiSq, df)]
  res1[, "onesidedp":=pNorm/2]
  
  pvals<-transpose(res1[,.(onesidedp)])
  setnames(pvals,paste0("pval_",n1))
  pval=res1[W==fparam$wlr,onesidedp]
  y<-cbind(pval,pvals)

  return(y)
}

#######################################################################################
## Landmark analysis
#######################################################################################


