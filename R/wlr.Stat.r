#' Weighted Logrank Test Statistic
#'
#' Weighted logrank adopted from survMisc package
#'
#' The one-sided p-Values from the following tests can be obtained based on input from \code{fparam$stdset} and \code{fparam$APPLE}
#' \describe{
#'  \item{1}{log-rank (part of standard set)} 
#'  \item{n}{Gehan-Breslow generalized Wilcoxon (part of standard set)}
#'  \item{sqrtN}{Tarone-Ware (part of standard set)}
#'  \item{S1}{Peto-Peto's modified survival estimate (part of standard set)}
#'  \item{S2}{modified Peto-Peto (by Andersen) (part of standard set)}
#'  \item{APPLE}{Simplified APPLE method with user specified late separation time}
#'  \item{FH(rho,gamma)}{Fleming-Harrington with user specified pairs of \code{rho} and \code{gamma}}
#'  }
#'
#' @param survival time-to-event variable
#' @param cnsr censoring variable: 1=censoring, 0=event
#' @param trt  treatment varaible. Accepted values are either "experimental" or "control"
#' @param stra stratification variable for stratified weighted logrank test. Default is \code{NULL}.
#' @param fparam a list input to request additional FH test and time of delayed separation in the simplified APPLE method
#' \describe{
#'  \item{\code{fparam$stdset}}{if not \code{NULL}, the standard set of tests are added} 
#'  \item{\code{fparam$rho}}{a vector of \code{rho} to be used in Fleming-Harrington FH(rho, gamma)} 
#'  \item{\code{fparam$gamma}}{a vector of \code{gamma} to be used in Fleming-Harrington \code{FH(rho, gamma)}. 
#'  The length of \code{gamma} needs to be the same as \code{rho}} 
#'  \item{\code{fparam$wlr}}{test name of which the p-Value will be stored in the variable \code{pval}. 
#'  Accepted values are "1", "n", "sqrtN", "S1", "S2", "APPLE", or "FH(rho,gamma)" where \code{rho} and \code{gamma} is one of the
#'  pair specified by user in the \code{fparam}}
#'  \item{\code{fparam$APPLE}}{If not left \code{NULL}, it is the time of separation used in the simplified piecewise exponential method.}
#'  } 
#' @return The function return a list with the follow components
#' \describe{
#'  \item{pval}{One-sided p-Value from user specified test} 
#'  \item{pval_[XXX]}{One-sided test from various weighted tests}
#'  }
#'   
#' @examples
#' # weighted logrank on the simulated data
#' library(survMisc)
#' medC = 6 
#' hr <- c(1, 0.6)
#' intervals <- 3 
#' gamma <- c(2.5, 5,  7.5,  10) ## a ramp-up enrollment
#' R     <- c(2  , 2,  2  ,  6 ) ## enrollment period: total of 12 months
#' eta <- -log(0.99) ## 1% monthly dropout rate
#' sim1 <- nphsim(nsim=2,lambdaC=log(2)/medC,lambdaE=log(2)/medC*hr, ssC=300,ssE=300,
#'                intervals=intervals,gamma=gamma, R=R,eta=eta)
#' test1 <- simtest(x=sim1, anaD=c(250,300), method=wlr.Stat,fparam=list(rho=c(0,0,1), 
#'                  gamma=c(0,1,1), wlr='FH(0,1)', APPLE=3))
#' test1$result[]
#' 
#' # direct function call (without cutoff)
#' wlr.Stat(surv=sim1$simd$survival, cnsr=sim1$simd$cnsr, trt=sim1$simd$treatment,
#'          fparam=list(rho=c(0,0,1), gamma=c(0,1,1), wlr='FH(0,1)', APPLE=3))
#' 
#' @export
#' @import data.table survival
#' @importFrom survMisc ten COV sf
#' 
wlr.Stat <- function(survival, cnsr, trt, stra = NULL, fparam) {
  if (length(fparam$rho)!=length(fparam$gamma)){
    stop("rho and gamma have to be of the same length")
  }
  if (is.null(stra)){stra=1}
  stralvl <- unique(stra)
  r<-NULL   # to hold the final statistics for each strata 
  for (i in 1:length(stralvl)) {
    d <- data.table(survival = survival, cnsr = cnsr, trt = trt, stra=stra)
    d <- d[stra==stralvl[i]]
    x <- ten(survfit(Surv(survival, 1 - cnsr) ~ trt, data = d))
    if (!attr(x, "sorted")=="t") setkey(x, t)
    t1 <- x[e > 0, t, by = t][, t]
   
    FHn <- paste("FH(", fparam$rho, ",", fparam$gamma,")", sep="")
   
    n1 <- c(FHn)
    
    if (!is.null(fparam$stdset)){
      n1 <- c(n1, c("1", "n", "sqrtN", "S1", "S2"))
    }
    
    if (!is.null(fparam$APPLE)){
      n1 <-c(n1, "APPLE")
    }
    wt <- data.table(array(data = 1, dim = c(length(t1), length(n1))))

    if (!is.null(fparam$wlr)){
      if (!fparam$wlr %in% n1){
        stop("fparam$wlr value is invalid. Refer to the help document for a list of allowed values.")
      }
    }
    
    data.table::setnames(wt, n1)
    
    
    if (!is.null(fparam$stdset)){
    ## Gehan-Breslow generalized Wilcoxon, weight = n
      data.table::set(wt, j = "n", value = x[e > 0, max(n), by = t][, V1])
      ## Tarone-Ware, weight = sqrt(n)
      data.table::set(wt, j = "sqrtN", value = wt[, sqrt(.SD), .SDcols = "n"])
      ## Peto-Peto, weight = S(t) = modified estimator of survival function
      data.table::set(wt, j = "S1", value = cumprod(x[e > 0, 1 - sum(e)/(max(n) + 1), by = t][, V1]))
      ## modified Peto-Peto (by Andersen), weight = S(t)n / n+1
      data.table::set(wt, j = "S2", value = wt[, S1] * x[e > 0, max(n)/(max(n) + 1), by = t][, V1])
    }
    
    ## Fleming-Harrington
    S3 <- sf(x = x[e > 0, sum(e), by = t][, V1], n = x[e > 0, max(n), by = t][, V1], what = "S")
    S3 <- c(1, S3[seq.int(length(S3) - 1L)])
    wt[, (FHn) := mapply(function(p, q) S3^p * ((1 - S3)^q), fparam$rho, fparam$gamma, SIMPLIFY=FALSE)]
  
    if (!is.null(fparam$APPLE)){
      ## simplified APPLE, w1=0 and w2=1, separate at delayed separation point
      data.table::set(wt, j = "APPLE", value = x[e > 0, ifelse(t < fparam$APPLE, 0, 1), by = t][, V1])
    }
    
    n2 <- c("W", "Q", "Var", "Z", "pNorm", "chiSq", "df", "pChisq")
    res <- data.table(matrix(0, nrow = ncol(wt), ncol = length(n2)))
    setnames(res, n2)
    set(res, j = 1L, value = n1)
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
    eMP1 <- unlist(eMP1[, .SD, .SDcols = (length(eMP1) - 1L)])
    data.table::set(res, j = "Q", value = colSums(wt * eMP1))
    data.table::set(res, j = "Var", value = colSums(wt^2 * cov1))
    r<-rbind(r,res)
  }  
  res<- r[,.(Q=sum(Q),Var=sum(Var)),by=W]
  res[, `:=`("Z", Q/sqrt(Var))]
  res[, `:=`("pNorm", 2 * (1 - stats::pnorm(abs(Z))))]
  res[, `:=`("chiSq", Q^2/Var)]
  res[, `:=`("df", 1)]
  res[, `:=`("pChisq", 1 - stats::pchisq(chiSq, df))]
  res[, `:=`("onesidedp", ifelse(Z < 0, pNorm/2, (1 - pNorm/2)))]
  
  pvals <- transpose(res[, .(onesidedp)])
  setnames(pvals, paste0("pval_", n1))
  if (!is.null(fparam$wlr)){
    pval = res[W == fparam$wlr, onesidedp]
    pvals <-cbind(pval, pvals)
  }
  y <- cbind(pvals)
  
  return(y)
}

