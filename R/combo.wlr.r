#' Weighted Logrank Combo Test and Estimation
#'
#' Description
#'
#' Details are provided here
#'
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment varaible. Accepted values are either "experiment" or "control"
#' @param stra stratification variable. Currently is set to \code{NULL}.
#' @param fparam a list input for additional function call arguments
#' \describe{
#'  \item{\code{fparam$rgs}}{required: a list of rho gamma pair in the form of \code{list(c(rho1, gamma1), c(rho2, gamma2), etc...)}} 
#'  \item{\code{fparam$draw}}{Number of resampling} 
#'  } 
#' @return The function return a list with the follow components
#' \describe{
#'  \item{pval}{One-sided p-Value corresponding to Zmax} 
#'  \item{rho, gamma}{rho and gamma corresponding to Zmax}
#'  \item{Zmax}{Zmax}
#'  }
#'   
#' @examples
#' medC = 6 
#' hr <- c(1, 0.6)
#' intervals <- 3 
#' gamma <- c(2.5, 5,  7.5,  10) ## a ramp-up enrollment
#' R     <- c(2  , 2,  2  ,  6 ) ## enrollment period: total of 12 months
#' eta <- -log(0.99) ## 1% monthly dropout rate
#' sim1 <- nphsim(nsim=2,lambdaC=log(2)/medC,lambdaE=log(2)/medC*hr, ssC=300,ssE=300,
#'                intervals=intervals,gamma=gamma, R=R,eta=eta)
#'                
#' rgs <- list(c(0, 0), c(0, 1), c(0, 2), c(0, 3), c(1, 0), c(2, 0), c(2, 2))    
#' draws <- 10        
#' test1<-simtest(x=sim1,anaD=c(250),method=combo.wlr,fparam=list(rgs=rgs,draws=draws))
#' test1$result[]
#' 
#' # direct function call (without cutoff)
#' combo.wlr(survival=sim1$simd$survival, cnsr=sim1$simd$cnsr, trt=sim1$simd$treatment,
#'           fparam=list(rgs=rgs,draws=draws))
#' 
#' @export
#' @import data.table survival
#' 


combo.wlr<- function(survival, cnsr, trt, stra = NULL, fparam) {
  temp <- fit.combo.wlr(rgs = fparam$rgs, time = survival, delta = 1-cnsr, z = (trt=="experimental"), draws = fparam$draw,
                        plot.rg = FALSE, print.results = FALSE, outfile = NULL)
  return(list(rho =   temp$CoxMax[1],
              gamma = temp$CoxMax[2],
              Zmax =  temp$CoxMax[3],
              pval =  temp$CoxMax[4],
              hr =    temp$CoxMax[5],
              hrL =   temp$CoxMax[6],
              hrU =   temp$CoxMax[7],
              hrL.bc =temp$CoxMax[8],
              hrU.bc =temp$CoxMax[9]
              ))
}

