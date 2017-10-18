#' Weighted Kaplan-Meier Test Statistic
#'
#' Produce the Weighted Kaplan-Meier test statistic. 
#' Test statistic is constructed based on the paper of Pepe and Fleming (1989).
#' Coding is based on work from a non-proportional hazards working group at Merck; original author ??
#' 
#' Although it is of interest that the statistic will reduce to RMST when the 
#' weight function = 1, the stability conditions proposed by Pepe and Fleming (1989) 
#' rule out a constant weight function when censoring is present. 
#' Input: dataframe 'indata' contains the following variables:
#' Event.T -> Event time (failure/censor); Status -> 0 = Censored, 1 = Event of Interest; 
#' Z -> 0 = Control, 1 = Active
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment varaible. Accepted values are either "experiment" or "control"
#' @param stra stratification variable. Default is \code{NULL} (currently not implemented)
#' @param fparam parameter description. Set to \code{NULL} (currently not used)
#' @return Code returns the Z statistic and p-value.
#' \describe{
#'  \item{pval}{One-sided p-Value from weighted Kaplan-Meier test} 
#'  \item{z}{test statistics}
#'  \item{Df}{degree of freedom. Currently always set to \code{NA}}
#'  }
#' @examples
#' # weighted Kaplan-Meier test on the simulated data
#' library(survival)
#' medC = 6 
#' hr <- c(1, 0.6)
#' intervals <- 3 
#' gamma <- c(2.5, 5,  7.5,  10) ## a ramp-up enrollment
#' R     <- c(2  , 2,  2  ,  6 ) ## enrollment period: total of 12 months
#' eta <- -log(0.99) ## 1% monthly dropout rate
#' sim1 <- nphsim(nsim=1,lambdaC=log(2)/medC,lambdaE=log(2)/medC*hr, ssC=300,ssE=300,
#'                intervals=intervals,gamma=gamma, R=R,eta=eta)
#' test1 <- simtest(x=sim1, anaD=c(250,300), method=wkm.Stat)
#' test1$result[]
#' 
#' # direct function call (without cutoff)
#' wkm.Stat(surv=sim1$simd$survival, cnsr=sim1$simd$cnsr, trt=sim1$simd$treatment)
#' @export
#' @import survival
wkm.Stat <- function(survival,cnsr,trt,stra=NULL,fparam=NULL){
  indata<-data.table(Event.T=survival,Status=1-cnsr,Z=(trt=='experimental'))
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
    pval = 1-pnorm(stat),
    z = stat,
    Df = NA
  ))
}