#' Update gsDesign object boundary
#' 
#' Update gsDesign object boundary based on actual number of events. 
#' Returns Z-values for upper or lower bounds for an updated design.
#' Internal function only.
#' 
#' Details to be added
#' 
#' @param dth Vector with number of events
#' @param d object of class 'gsDesign'
#' @param b 1 for upper bound, 0 for lower bound
#' 
#' @keywords internal
gsUpdate <- function(dth, d, b) {
  ##b=1 upper bound, b=0 lower bound
  dnew <-
    gsDesign(
      k = d$k,
      test.type = d$test.type,
      alpha = d$alpha,
      beta = d$beta,
      n.fix = d$n.fix,
      sfu = d$upper$sf,
      sfupar = d$upper$param,
      sfl = d$lower$sf,
      sflpar = d$lower$param,
      n.I = dth,
      maxn.IPlan = d$n.I[d$k]
    )
  if (b == 1) {
    y <- pnorm(dnew$upper$bound, lower.tail = FALSE)
  } else{
    y <- pnorm(dnew$lower$bound, lower.tail = FALSE)
  }
  return(y)
}


