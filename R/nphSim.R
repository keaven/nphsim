#' Simulate a Clinical Trial with Piecewise Exponential Time-to-Event Outcomes
#'
#' Simulate two-arm time-to-event data using the piecewise exponential distribution \code{rpwexp()}. 
#' User can specify enrollment speed as well as drop out rate separately for each arm. Additionaly if user has created
#' a gsSurv object from \code{\link[gsDesign:nSurv]{gsDesign}} it can be used as input to supply simulation parameters. The only
#' censoring mechanism is from dropout of the study and no administrative censoring is implemented.
#'
#' All the simulation parameters: sample size, hazard rate in each interval and the interval duration, 
#' enrollment time period and enrollment speed, and the piecewise dropout rate for the same interval duration need to be provided 
#' unless a gsSurv object is provided, in which case the individual parameter can be left as \code{NULL}. If not \code{NULL} the value will
#' overwrite the corresponding value from the gsSurv object.  
#'
#' @param nsim Number of simulations
#' @param lambdaC Hazard rate of control arm. Specify a vector for piecewise hazard with duration specified in "intervals"
#' @param lambdaE Hazard rate of experiment arm. Specify a vector for piecewise hazard with duration specified in "intervals"
#' @param intervals Duration of period in which hazard is constant as specified in lambdaC. A vector with length(lambdaC)-1
#' @param ssC Sample size of control arm
#' @param ssE Sample size of experiment arm
#' @param gamma A vector of rate of enrollment per unit of time
#' @param R A vector of duration of time periods for recruitment with rates specified in gamma; should be same length as gamma
#' @param fixEnrollTime if \code{TRUE} the enrollment period \code{R} is fixed and enrollment rate is adjusted proportionally to meet the sample size; 
#'  otherwise the last interval of R is adjusted to meet the sample size without adjustment to the enrollment rate.
#' @param eta A vector for dropout rate per unit time for control arm
#' @param etaE A vector for dropout rate per unit time for experiment arm
#' @param d A gsSurv object as input. The other inputs overwrite the corresponding parameters if not NULL
#' 
#' @return The function return a list with the follow components
#' \describe{
#'  \item{nsim, lambdaC, lambdaE, ssC, ssE, intervals}{as Input} 
#'  \item{gamma}{actual enrollment rate per enrollment period. If \code{fixEnrollTime} is not \code{TRUE} this is the same as the input \code{gamma}}
#'  \item{R}{actual enrollment period. If \code{fixEnrollTime} is \code{TRUE} this is the same as input \code{R}, otherwise
#'  the last interval of the input \code{R} will be adjusted to meet the specified sample size while keeping \code{gamma} fixed.}
#'  \item{etaC}{same as eta} 
#'  \item{etaE}{as Input or equal to etaC if Input is NULL} 
#'  \item{simd}{data table object that stores the simulated data
#'    \itemize{
#'    \item{sim: simulation sequence number}
#'    \item{treatment: "control" or "experiment"}
#'    \item{enterT: calendar time a subject enters the study}
#'    \item{ct: calendar time of event/censoring. Equals to \code{enterT + survival}}
#'    \item{survival: simulated time-to-event value}
#'    \item{cnsr: censoring status. 1 = censored, 0 = event}
#'    }
#'  } 
#'  }
#' 
#' 
#' @examples
#' # Simulate a two-arm study with overall survival endpint using proportional harzard of 0.7 and median survival of 6 months. 
#' # The number of subjects in each arm is 300.
#' library(survival)
#' medC = 6 
#' hr <- 0.7
#' intervals <- NULL 
#' gamma <- c(2.5, 5,  7.5,  10) ## a ramp-up enrollment
#' R     <- c(2  , 2,  2  ,  6 ) ## enrollment period: total of 12 months
#' eta <- -log(0.99) ## 1% monthly dropout rate

#' sim1<-nphsim(nsim=10,lambdaC=log(2)/medC,lambdaE=log(2)/medC*hr, ssC=300,ssE=300,intervals=intervals,gamma=gamma, R=R,eta=eta)
#' km<-survfit(Surv(survival,1-cnsr)~treatment,data=sim1$simd[sim==1])
#' plot(km,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36),mark.time=TRUE)
#' ll <- gsub("x=","",names(km$strata))  ## legend labels
#' legend("top",legend=ll,lty=1:2,horiz=FALSE,bty='n')
#'  
#' 
#' # Simulate survival endpoint with delayed separation of KM curve: HR=1 for the first 3 months and 0.65 afterwards 
#' hr <- c(1,0.65)
#' intervals <- 3 
#' sim2<-nphsim(nsim=10,lambdaC=log(2)/medC,lambdaE=log(2)/medC*hr, ssC=300,ssE=300,intervals=intervals,gamma=gamma, R=R,eta=eta)
#' km<-survfit(Surv(survival,1-cnsr)~treatment,data=sim1$simd[sim<=10])
#' plot(km,xlab="Month Since Randomization",ylab="Survival",lty=1:2,xlim=c(0,36))
#' ll <- gsub("x=","",names(km$strata))  ## legend labels
#' legend("top",legend=ll,lty=1:2,horiz=FALSE,bty='n')
#' 
#' 
#' # Use gsSurv as input
#' library(gsDesign)
#' gs <- gsSurv ( k = 3, test.type = 4, alpha = 0.025, beta = 0.05, timing = c( 0.5,0.75 ), sfu = sfHSD , sfupar = c( -4 ), sfl = sfHSD, sflpar = c( -12 ), lambdaC = log(2) / 6, hr = 0.65, hr0 = 1, eta = 0.01, gamma = c( 2.5,5,7.5,10 ), R = c( 2,2,2,6 ) , S = NULL , T = 15 , minfup = 3 , ratio = 1) 
#' sim3 <- nphsim(nsim=10,d=gs) 

#' @export
#' @import gsDesign data.table
nphsim <- function(nsim = 100
                  ,lambdaC = NULL
                  ,lambdaE = NULL
                  ,intervals = NULL
                  ,ssC = NULL
                  ,ssE = NULL
                  ,gamma = NULL
                  ,R = NULL
                  ,fixEnrollTime = TRUE
                  ,eta = NULL
                  ,etaE = NULL
                  ,d = NULL
                  )
{
  ## use gsSurv object if provided, values can be overwritten by other parameters
  if (!is.null(d)) {
    if (!is(d, "gsDesign")) {
      stop("d should be an object of class gsDesign")
    }
    ssC <- if(is.null(ssC)){ceiling(max(d$eNC))}else {ssC}
    ssE <- if(is.null(ssE)){ceiling(max(d$eNE))}else {ssE}
    lambdaC <- if(is.null(lambdaC)){d$lambdaC}else{lambdaC}
    lambdaE <- if(is.null(lambdaE)){d$lambdaC*d$hr}else{lambdaE}
    if (length(lambdaC)==1) {lambdaC=rep(lambdaC,length(lambdaE))}
    if (length(lambdaC)==1){
      intervals <- NULL
    } else {
      if (is.null(intervals)){
        intervals <- d$S
      }
    }
    gamma <- if(is.null(gamma)){d$gamma}else{gamma}
    R <- if(is.null(R)){d$R}else{R}
    eta <- if(is.null(eta)){d$etaC}else{eta}

  } else {  # if gsSurv object is not provided, none of the other parameters can be NULL
    if (is.null(lambdaC) | is.null(lambdaE) | is.null(ssC) | is.null(ssE) | is.null(gamma)
        | is.null(R) | is.null(eta)){
      stop("All arguments need to be provided if gsSurv object is not specified in d.")
    }
  }
  tnum <- (ssC+ssE)*nsim
  etaC <- eta
  etaE <- if(is.null(etaE)){etaE=etaC}else{etaE}

  ## simulate survival based on piece-wise hazard rate
  tC <- rpwexp(ssC*nsim, rate = lambdaC, intervals = intervals)
  tE <- rpwexp(ssE*nsim, rate = lambdaE, intervals = intervals)
  ## adding LTFU
  xC<-data.table(sim=c(1:nsim),t=tC,treatment='control')
  xC[,ltfuT:=rpwexp(ssC*nsim, rate = etaC+1e-8, intervals = intervals)]
  xE<-data.table(sim=c(1:nsim),t=tE,treatment='experiment')
  xE[,ltfuT:=rpwexp(ssE*nsim, rate = etaE+1e-8, intervals = intervals)]

  x <- rbind(xC,xE)
  x[,t1:=ifelse(t>ltfuT, ltfuT, t)]
  x[,cnsr1:=ifelse(t>ltfuT, 1, 0)]

  ## uniform enrollment in each intervals of R
  if (isTRUE(fixEnrollTime)){ # adjust enrollment rate while keeping R fixed
    aR <- R # keep actual enrollment period the same as input R
    agamma <- (ssC+ssE)*gamma/sum(gamma*R)  # Actual enrollment rate
    x[,enterT:=sample(rpwexp(.N,rate=.N*gamma/sum(gamma*R),intervals=R[1:length(gamma)-1],cumulative=TRUE)),by=sim]
  } else{ # adjust the last interval of R while keeping enrollment rate fixed
    agamma <- gamma
    aR <- R
    if ((ssC+ssE)-sum(gamma[1:length(gamma)-1]*R[1:length(gamma)-1])<=0) {
      stop ("User requested to keep the enrollment rate fixed but the total enrollment is already greater than the specified sample size prior to the last interval of R.")
    } else{
      aR[length(gamma)] <- ((ssC+ssE)-sum(gamma[1:length(gamma)-1]*R[1:length(gamma)-1]))/gamma[length(gamma)]
      x[,enterT:=sample(rpwexp(.N,rate=gamma,intervals=R[1:length(gamma)-1],cumulative=TRUE)),by=sim] 
    }
  }
  

  ## ct: calendar time a subject had event/censoring
  x[,ct:=t1 + enterT]
  x[,survival:=t1][,cnsr:=cnsr1]
  x[,c("cnsr1","t1","t","ltfuT"):=NULL]  ## remove intermediate variables
  x[,treatment:=relevel(factor(treatment),ref='control')]  ## set control as reference level

  y<-list(nsim=nsim,
          lambdaC=lambdaC,
          lambdaE=lambdaE,
          ssC=ssC,
          ssE=ssE,
          intervals=intervals,
          gamma=agamma,
          R=aR,
          etaC=etaC,
          etaE=etaE,
          simd=x)
  return(y)
}
#'
#' @importFrom(gsDesign)




