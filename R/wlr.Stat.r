#' Weighted Logrank Test Statistic
#'
#' Weighted logrank adopted from survMisc package
#'
#' Details text
#'
#' @param survival parameter description
#' @param cnsr parameter description
#' @param trt  parameter description
#' @param stra parameter description
#' @param fparam parameter description
#' @examples
#' # TBD
#' @export
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
#'
#' @import(doParallel)
#' @import(plyr)
#' @import(reshape2)
#' @import(survival)
#' @import(survMisc)
