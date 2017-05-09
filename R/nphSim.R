nphsim <- function(nsim = 100 # Number of simulations
                  ,lambdaC = NULL # Hazard rate of control arm. Specify a vector for piecewise hazard with duration specified in "intervals"
                  ,lambdaE = NULL # Hazard rate of experiment arm. Specify a vector for piecewise hazard with duration specified in "intervals"
                  ,intervals = NULL # Duration of period in which hazard is constant as specified in lambdaC. A vector with length(lambdaC)-1
                  ,ssC = NULL # Sample size of control arm
                  ,ssE = NULL # Sample size of experiment arm
                  ,gamma = NULL # A vector of rate of enrollment in unit time
                  ,R = NULL # A vector of duration of time periods for recruitment with rate specified in gamma
                  ,eta = NULL # A scalar for dropout rate per unit time for control arm
                  ,etaE = NULL # A scalar for dropout rate per unit time for experiment arm
                  ,d = NULL # An gsSurv object as input. The other inputs overwrite the corresponding parameters if not NULL
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
    if (is.null(lambdaC) | is.null(lambdaE) | is.null(ssC) | is.null(ssE) | is.null(gamma) | is.null(R) | is.null(eta) | is.null(T)){
      stop("All arguments need to be provided if gsSurv objects is not specified in d.")
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
  x[,enterT:=sample(rpwexp(.N,rate=.N*gamma/sum(gamma*R),intervals=R[1:length(R)-1],cumulative=TRUE)),by=sim]

  ## ct: calendar time a subject had event/censoring
  x[,ct:=t1 + enterT]
  ## administrative censor at T
  ## keep this checking for now and will remove eventually
  T <- 99999
  x[,survival:= ifelse(ct>T, T-enterT, t1)]
  x[,cnsr:=ifelse(ct>T, 1, cnsr1)]
  x[,c("cnsr1","t1","t","ltfuT","ct"):=NULL]  ## remove intermediate variables
  x[,treatment:=relevel(factor(treatment),ref='control')]  ## set control as reference level

  y<-list(nsim=nsim,
          lambdaC=lambdaC,
          lambdaE=lambdaE,
          ssC=ssC,
          ssE=ssE,
          intervals=intervals,
          gamma=gamma,
          R=R,
          etaC=etaC,
          etaE=etaE,
          simd=x)
  return(y)
}





