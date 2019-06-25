#' Testing for a time-to-event study simulation
#' 
#' Generate additional simulated data at user specified analysis sequences (either driven by calendar time or number of events). 
#' Create simulation level summary of analysis timing and number of events, as well as testing results and boundary crossing 
#' 
#' \code{simtest()} takes the returned list from \code{nphsim()} and creates and analyzes the simulated data at various analysis time points 
#' depending on user input:
#' \itemize{
#'   \item{anatype = 'calendar': by specified calendar time with anaT}
#'   \item{anatype = 'event': by specified number of events with anaD}
#'   \item{anatype = 'both': use the later analysis time points of either anaT or anaD}
#' }
#' Different statistical testing and summarizing procedures can be specified in the method parameter. 
#' User defined statistical procedures can be used so long as it follows certain order of parameter setup in function definition.
#' Similar to \code{nphsim()}, \code{simtest()} can also take a gsDesign object as input, same or different from the one used in nphsim. 
#' If used, boundary crossing is checked against the boundaries in the gsDesign object at each analysis (either specified time points or events). 
#' The boundaries are re-calculated if the number of events at an analysis is different from the design.
#'
#' @param x list genearted by simulation
#' @param anaT A vector or matrix of calendar time for reporting, testing, and boundary crossing 
#' probability when a gsDesign object is specified.
#' @param anaD A vector of number of events for reporting, testing, and boundary crossing 
#' probability when a gsDesign object is specified.
#' @param anatype use specified timing ('anaT') if 'calendar'; 
#' use specified number of events ('anaD') if 'event'; 
#' use the maximum of the two if 'both'
#' @param method statistical testing method to be applied to simulated data. 
#' @param stratum the variable name in the simulated data for the stratified analysis 
#' @param fparam additional parameters needed in the user defined testing functions. This has to match the parameters in the testing function.
#' @param d A 'gsDesign' object. Timing/Events and the boundary will be updated upon output.
#'
#' @export
#' @return The function return a list with the follow components
#' \describe{
#'  \item{anaT, anaD, anatype, method}{as Input} 
#'  \item{gsobj}{same as d, the gsDesign input} 
#'  \item{simd}{data table object that stores the simulated data by analysis sequence}
#'    \itemize{
#'    \item{sim: simulation sequence number}
#'    \item{analysis: analysis sequence corresponding to first interim, second interim, etc...}
#'    \item{t: calendar time corresponding to each analysis sequence}
#'    \item{D: number of events corresponding to each analysis sequence}
#'    \item{treatment: "control" or "experiment"}
#'    \item{enterT: calendar time a subject enters the study}
#'    \item{aval: simulated time-to-event value}
#'    \item{cnsr: censoring status. 1 = censored, 0 = event}
#'    }
#'  \item{result}{data table object that stores the simulation level summaries}
#'    \itemize{
#'    \item{sim: simulation sequence number}
#'    \item{analysis: analysis sequence corresponding to first interim, second interim, etc...}
#'    \item{t: calendar time corresponding to each analysis sequence}
#'    \item{D: number of events corresponding to each analysis sequence}
#'    \item{NE, NC: Number of subjects for experiment/control group}
#'    \item{DE, DC: Number of events for experiment/control group}
#'    \item{pval: p-Value from testing if \code{method} parameter is specified}
#'    \item{pvupper, pvlower: upper/lower bound nominal p-Value from the group sequential design if \code{method} and \code{d} parameters are specified)}
#'    \item{xeff, xfut: upper/lower bound crossing status: 1=Yes,0=No, if \code{method} and \code{d} parameters are specified). 
#'    xeff and xfut take into account the previous boundary crossing status for each simulation. 
#'    E.g. if upper boundary was acrossed in a previous analysis, all subsequent xeff will be set to 1. 
#'    Similarly if lower boundary was acrossed in a previous analysis, all subsequent xeff will be set to 0.}
#'    \item{addtional values returned by the testing functions}
#'    }
#'  }
#'  
#' @examples
#' # Use a gsSurv object as both input for simulation and testing. 
#' # A logrank p-value and HR from cox model is reported.
#' library(gsDesign)
#' gs <- gsSurv (k = 3, test.type = 4, alpha = 0.025, beta = 0.05, timing = c( 0.5,0.75 ), 
#'               sfu = sfHSD , sfupar = c( -4 ), sfl = sfHSD, sflpar = c( -12 ), 
#'              lambdaC = log(2) / 6, hr = 0.65, hr0 = 1, eta = 0.01, 
#'              gamma = c( 2.5,5,7.5,10 ), R = c( 2,2,2,6 ) , S = NULL , T = 15 , minfup = 3 , ratio = 1) 
#' sim1 <- nphsim(nsim=10,d=gs)
#' test1 <- simtest(x=sim1,anatype='event',method='LR', d=gs)
#' test1$result
#' plotsim(test1$result,y=c("hr"),dg=2,yt="Hazard Ratio",b=1,v=c(2,0.8))
#' 
#' 
#' # specify fixed analysis time
#' test2 <- simtest(x=sim1,anatype='calendar',anaT=c(10, 13, 16), method='LR', d=gs)
#' plotsim(test2$result,y=c("hr"),dg=2,yt="Hazard Ratio",b=1,v=c(2,0.8))
#' 
#' # specify by-simulation analysis time
#' anaT<-matrix(c(9,11,15,10,12,16),byrow=TRUE,nrow=00,ncol=3)
#' test3 <- simtest(x=sim1,anatype='calendar',anaT=anaT, method='LR', d=gs)
#' 
#' 
#' # specify number of events
#' test4 <- simtest(x=sim1,anatype='event',anaD=c(140, 200, 290), method='LR', d=gs)
#'
#' # A user defined testing function. Additional information needed for the testing function 
#' # can be included in the fparam. 
#' # and provided in the \code{simtest() fparam} parameter.
#' library(survRM2)
#' test.usr1<-function (survival,cnsr,trt,stra=NULL,fparam=NULL) {
#' a<-rmst2(time=survival,status=1-cnsr,arm=(trt=='experiment'),tau=fparam)
#' b<-a$unadjusted.result
#' pval <- b[1,4]
#' y<-list(pval=round(pval,4), est=b[1,1], estlb=b[1,2],estub=b[1,3])
#' return(y)
#' }
#' test5 <- simtest(x=sim1,anatype='event',anaD=c(200, 290), method=test.usr1, fparam=9)
#' plotsim(test5$result,y=c("estlb","estub"),dg=2,yt="Estimate",v=c(2,0.8))
#' 
#' @import gsDesign data.table
#######################################################################################
simtest <- function(x 
                    ,anaT=NULL 
                    ,anaD=NULL 
                    ,anatype='event'
                    ,method=NULL
                    ,stratum=NULL
                    ,fparam=NULL
                    ,d=NULL 
                    )
{
  if (anatype %in% c('event', 'calendar', 'both')) {
  } else{
    stop("anatype needs to be 'calendar', 'event', or 'both'.")
  }  
  ## use gsSurv object if provided
  flag <- NULL
  if (!is.null(d)) {
    if (!is(d, "gsDesign")) {
      stop("d should be an object of class gsDesign")
    }
    if (is.null(anaD) &
        is.null(anaT)) { ## Use the original gsDesign boundary
    } else{
      flag = 1 ## gsDesign boundary needs to be updated when flag=1
    }
    anaT <- if (is.null(anaT)) {
      d$T
    } else{
      anaT
    }
    anaD <- if (is.null(anaD)) {
      c(ceiling(d$eDC + d$eDE))
    } else{
      anaD
    }
  }

  if (anatype=='calendar'){
    ## use anaT for calendar based interim analyses
    sim=c(1:x$nsim)
    if (is.vector(anaT)){ ## same analysis time for all simulations 
      anaTDT<-data.table(sim,matrix(anaT,ncol=length(anaT)))
    }
    if (is.matrix(anaT)){ ## vary analysis time for each simulation 
      anaTDT<-data.table(sim,anaT)
    }
    DT.t<-melt(anaTDT, id.vars=c("sim"))
    DT.t<-DT.t[,.(sim=sim,t=value,analysis=as.numeric(substr(variable,2,10)))]
    DT.sim <- merge(DT.t, x$simd, by = "sim", allow.cartesian = TRUE)
    DT.sim [, aval := ifelse(survival + enterT > t, t - enterT, survival)
            ][, cnsr := ifelse(survival + enterT > t, 1, cnsr)
              ][, c("survival"):=NULL]
  } else if (anatype == 'event') {
    ## use anaD for event driven interim analyses
    dt <- data.table(x$simd,
                     #ct = x$simd$survival + x$simd$enterT, # ct is added in nphsim returned value so this is not needed
                     k = 1)
    setorderv(dt, c("sim", "cnsr", "ct"), c(1, 1, 1))
    dt <- dt[, evn := rep(c(1:(x$ssC + x$ssE)),x$nsim)) # replaced by Ray Lin: c(1:(x$ssC + x$ssE))
             ][, evn := ifelse(cnsr == 1, 9999, evn)]
    DT.D <- data.table(D = anaD,
                       analysis = c(1:length(anaD)),
                       k = 1)
    DT.sim1 <- merge(dt, DT.D, by = "k", allow.cartesian = TRUE)
    survt <- DT.sim1[evn == D, .(sim, t = ct, D)]
    DT.sim <- merge(DT.sim1, survt, by = c("sim", "D"), all.x = TRUE)
    DT.sim[, t := ifelse(is.na(t), 99999, t)
           ][, survival2 := ifelse(ct <= t, survival, t - enterT)
             ][, cnsr2 := ifelse(ct <= t, cnsr, 1)
               ][, aval := survival2][, cnsr := cnsr2
                                      ][, c("k","survival","evn","survival2","cnsr2","ct"):=NULL]
    
  } else if (anatype == 'both') {## take the max event/calendar time of anaT and anaD for each analysis
    if (!is.vector(anaT)){
      stop("When anatype is set to 'both', anaT has to be a vector.")
    }
    if (length(anaT)!=length(anaD)){
      stop("anaT and anaD need to have the same length")
    } else {
      dt <- data.table(x$simd,
                       k = 1)
      setorderv(dt, c("sim", "cnsr", "ct"), c(1, 1, 1))
      dt <- dt[, evn := c(1:(x$ssC + x$ssE))
               ][, evn := ifelse(cnsr == 1, 9999, evn)]
      DT.D <- data.table(D = anaD,
                         analysis = c(1:length(anaD)),
                         k = 1)
      DT.sim1 <- merge(dt, DT.D, by = "k", allow.cartesian = TRUE)
      survt <- DT.sim1[evn == D, .(sim, anaDt = ct, D)]
      DT.sim <- merge(DT.sim1, survt, by = c("sim", "D"), all.x = TRUE)
      DT.sim[, anaDt := ifelse(is.na(anaDt), 99999, anaDt)]
      DT.sim[, anaTt:=anaT[analysis]][,t:=max(anaDt,anaTt),by=.(sim, analysis)]
      DT.sim[, survival2 := ifelse(ct <= t, survival, t - enterT)
             ][, cnsr2 := ifelse(ct <= t, cnsr, 1)
               ][, aval := survival2][, cnsr := cnsr2
                                      ][, c("k","survival","evn","survival2","cnsr2","ct"):=NULL]
    }
  }
  
  DT.sim<-DT.sim[aval > 0] ## remove subjects enrolled after specified time or number of events
  setorderv(DT.sim, c("sim", "analysis"), c(1, 1))
  tmp1 <- DT.sim[treatment == 'experimental', .(t = t[.N],NE = .N,DE = sum(1 - cnsr)), by = .(sim, analysis)]
  tmp2 <- DT.sim[treatment == 'control', .(NC = .N, DC = sum(1 - cnsr)), by = .(sim, analysis)]
  tmp <- merge(tmp1, tmp2, by = c("sim", "analysis"), all = TRUE)[, D :=DC + DE]
  
  result<-NULL
  if (!is.null(method)){  ## start testing procedures
    e=substitute(stratum) ## use user provided stratum variable 
    
    if (is.character(method)) {  ## if a pre-defined method is used
      if (!is.element(method,c("LR"))){
        stop("Character specification of method may only be LR")
      }
      result<-DT.sim[, c(test.lr(aval,cnsr,treatment,eval(e))), by=.(sim,analysis)]  
    } else if (is.function(method)){ ## if user defined function is used
      result<-DT.sim[, c(method(aval,cnsr,treatment,eval(e),fparam)), by=.(sim,analysis)]  
      if (!any(names(result)=="pval")){
        warning("user defined function did not return a p-Value with name 'pval'. The boundary crossing cannot be checked.")
      }
    }
    
    result<-merge(tmp,result,by=c("sim","analysis"))
    if (!is.null(d)) {
      if (is.null(flag)){  ### use the original boundary nominal p-values
        result[,pvupper:=pnorm(d$upper$bound,lower.tail=FALSE)
               ][,pvlower:=pnorm(d$lower$bound,lower.tail=(d$test.type==2))]
      } else {  ### update boundary nominal p-values
        result[,pvupper:=gsUpdate(D,d,1),by=sim
               ][,pvlower:=gsUpdate(D,d,2),by=sim]
      }
      ## checking boundary crossing while taking history into account
      result[,xeff1:=ifelse(pval<pvupper,1,0)
             ][,xfut1:=ifelse(pval>pvlower,1,0)
               ][,xeff2:=ifelse(cumsum(xeff1)>0,1,0), by=sim
                 ][,xfut2:=ifelse(cumsum(xfut1)>0,1,0), by=sim
                   ][,fstfl:=ifelse(sum(xeff2)>sum(xfut2),1,0), by=sim
                     ][,xeff:=xeff2*fstfl][,xfut:=xfut2*(1-fstfl)
                                           ][,c("fstfl","xeff1","xeff2","xfut1","xfut2"):=NULL]
    }
  }
  
  
  if (is.null(result)){resultall<-tmp}else{resultall<-result}
  
  y<-list(anaT=anaT, anaD=anaD, anatype=anatype, method=method,result=resultall, gsobj=d,simd=DT.sim)  
  
  
  return(y)
}  
