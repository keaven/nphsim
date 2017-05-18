#' Testing for a time-to-event study simulation
#' 
#' Description here
#' simulation test: provide # of events/timing at specified timing/# of events
#'                  provide HR and p-value at specified timing/events
#'                  provide boundary crossing probability (experimental: cannot update boundary at different # of events)
#' 
#' Details here
#'
#' @param x list genearted by simulation
#' @param anaT Calendar time for reporting, testing, and boundary crossing 
#' probability when a gsDesign object is specified.
#' @param anaD Events driven for reporting, testing, and boundary crossing 
#' probability when a gsDesign object is specified.
#' @param anatype use specified timing ('anaT') if 'calendar'; 
#' use specified number of events ('anaD') if 'event'; 
#' use the maximum of the two if 'both'
#' CAN INDIVIDUAL VALUES BE MISSING SO that, e.g. SOME CAN BE BOTH, BUT FINAL IS EVENTS?
#' @param method statistical testing method to be applied to simulated data. 
#' @param stratum definition TBD
#' @param fparam definition TBD
#' @param d=NULL A 'gsDesign' object. Timing/Events and the boundary will be 
#' updated upon output.
#' @examples
#' # examples here
#'
#' @export
#' @return
#' Needs update!
#' @import gsDesign
#' @import data.table
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
    DT.t <- data.table(t = anaT,
                       analysis = c(1:length(anaT)),
                       k = 1)
    DT.x <- data.table(x$simd, k = 1)
    DT.sim <- merge(DT.t, DT.x, by = "k", allow.cartesian = TRUE)
    DT.sim [, aval := ifelse(survival + enterT > t, t - enterT, survival)
            ][, cnsr := ifelse(survival + enterT > t, 1, cnsr)
              ][, c("k","survival"):=NULL]
  } else if (anatype == 'event') {
    ## use anaD for event driven interim analyses
    dt <- data.table(x$simd,
                     ct = x$simd$survival + x$simd$enterT,
                     k = 1)
    setorderv(dt, c("sim", "cnsr", "ct"), c(1, 1, 1))
    dt <- dt[, evn := c(1:(x$ssC + x$ssE))
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
    if (length(anaT)!=length(anaD)){
      stop("anaT and anaD need to have the same length")
    } else {
      dt <- data.table(x$simd,
                       ct = x$simd$survival + x$simd$enterT,
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
  tmp1 <- DT.sim[treatment == 'experiment', .(t = t[.N],NE = .N,DE = sum(1 - cnsr)), by = .(sim, analysis)]
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
        stop("user defined function must return a p-Value with name 'pval'")
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
