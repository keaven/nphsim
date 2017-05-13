#' Restricted Mean Survival Test Statistic
#'
#' RMST method with user specified trucation time.
#' Example of user provided testing function.
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
rmst.Stat<-function (survival,cnsr,trt,stra=NULL,fparam=NULL) {
  a<-rmst2(time=survival,status=1-cnsr,arm=(trt=='experiment'),tau=fparam$tt)
  b<-a$unadjusted.result
  pval <- b[1,4]/2
  y<-list(pval=pval, est=b[1,1], estlb=b[1,2],estub=b[1,3])
  return(y)
}
#' @import(survRM2)
