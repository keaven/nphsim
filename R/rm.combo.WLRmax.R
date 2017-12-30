#' Weighted Logrank combination with asymptotic p-value adjustment
#'
#' Description
#'
#' This routine provides several options for correcting p-values for the maximum of multiple 
#' weighted logrank test statistics.
#' The most favorable of these (i.e., minimally sufficient p-value correction) is provided by the "asymp" option which
#' uses a published asymptotic correlation between different weighted logrank tests (Karrison, 2016).
#' However, there are also more conservative (less powerful) options to correct using 
#' Holm, Hochberg, Hommel or Bonferroni corrections.
#'
#' @param time time-to-event variable
#' @param status event indicator: 0=censored, 1=event
#' @param arm treatment group indicator
#' @param wt a list with pairs of non-negative numbers setting rho and gamma for each 
#' Fleming-Harrington weighting sheme to be used
#' @param adjust.methods methods used to adjust Type I error for testing multiple weighted logrank tests
#' @param ties.method methods for dealing with tied event times
#' @param one.sided logical indicating 1-sided (TRUE) vs 2-sided (FALSE) testing; as currently written, the 
#' @param max indicator of extent of output (TRUE for maximum output)
#' @param alpha level for testing 
#' @return The function returns a list with the follow components
#' \describe{
#'  \item{pval}{One-sided p-Value corresponding to Zmax} 
#'  \item{rho, gamma}{rho and gamma corresponding to Zmax}
#'  \item{Zmax}{Zmax}
#'  }
#'
#' @details 
#' Note that the user needs to load the `dplyr` and `Matrix` packages since these have NOT been added to the imports list for the package.
#' The correlation between different weighted logrank test statistics noted in the Karrison (2016) is the same as the variance of a weighted
#' logrank statistic with the sum of rho values and sum of gamma values as the respective rho and gamma parameters.
#' @examples
#'
#' data(Ex1delayedEffect)
#' with(Ex1delayedEffect,rm.combo.WLRmax(time=month,status=evntd,arm=trt,adjust.methods="asymp",
#'                                    wt=list(a1=c(0,1),a2=c(1,1)))) 
#' with(Ex1delayedEffect,rm.combo.WLRmax(time=month,status=evntd,arm=trt,
#'                                    wt=list(a1=c(0,1),a2=c(1,1)))) 
#'                                    
#' @import survival survMisc
#'
#' @author Satrajit Roychoudhury, \email{satrajit.roychoudhury@@pfizer.com} 
#' @references TG Karrison (2016), Versatile tests for comparing survival curves based on weighted log-rank statistics, 
#' Stata Journal, 16:3 pp 678-690
#' @export
#' 
rm.combo.WLRmax<- function(time        = NULL,
                          status      = NULL,
                          arm         = NULL,
                          wt          = NULL,
                          adjust.methods = c("holm", "hochberg", "hommel", "bonferroni","asymp")[3],
                          ties.method = c("exact", "breslow", "efron")[2],
                          one.sided   = FALSE,
                          max         = TRUE,
                          alpha       = 0.025
)
{ 
  data.anal <- data.frame(cbind(time,status,arm))
  fit<- survMisc::ten(survival::survfit(survival::Surv(time, status) ~ arm, data = data.anal))  
  
  #Testing
  
  comp(fit, p= sapply(wt, function(x){x[1]}), q= sapply(wt, function(x){x[2]}))
  
  tst.rslt <- attr(fit, 'lrt')
  
  #Combination test (exact form)
  
  if(max & adjust.methods != "asymp"){
    tst.rslt1 <- rbind(tst.rslt[1,],subset(tst.rslt, grepl("FH", tst.rslt$W)))
    Z.tst.rslt1 <- tst.rslt1$Z
    if(one.sided){p.unadjusted <- pnorm(q=tst.rslt1$Z)}
    if(!one.sided){p.unadjusted <- 1- pnorm(q=abs(tst.rslt1$Z)) + pnorm(q=-abs(tst.rslt1$Z))}
    
    pval.adjusted <- p.adjust(p.unadjusted, method = adjust.methods)
    pval <- min(pval.adjusted)
    max.index <- which(p.unadjusted == min(p.unadjusted), arr.ind = TRUE)
  }
  
  
  if(max & adjust.methods == "asymp"){
    
    #Calculating the covariace matrix
    
    tst.rslt1 <- rbind(tst.rslt[1,],subset(tst.rslt, grepl("FH", tst.rslt$W)))
    
    Z.tst.rslt1 <- tst.rslt1$Z 
    q.tst.rslt1 <- tst.rslt1$Q 
    var.tst.rslt1 <- tst.rslt1$Var
    
    wt1 <- c(list(a0=c(0,0)), wt)
    combo.wt <- combn(wt1,2)
    
    combo.wt.list <- list()
    for(i in 1:ncol(combo.wt)){combo.wt.list[[i]] <- combo.wt[,i]}
    
    combo.wt.list.up <- lapply(combo.wt.list,function(a){mapply('+',a)})
    
    wt2 <- lapply(combo.wt.list.up, function(a){apply(a,1,'sum')/2})
    d1 <- data.frame(do.call(rbind,wt2))
    
    wt3 <- unique(wt2)
    d2 <- data.frame(do.call(rbind,wt3))
    
    fit2<- survMisc::ten(survival::survfit(survival::Surv(time, status) ~ arm, data = data.anal))  
    
    #Testing (for calculating the covariances)
    
    comp(fit2, p= sapply(wt3, function(x){x[1]}), q= sapply(wt3, function(x){x[2]}))
    
    tst.rsltt <- attr(fit2, 'lrt')
    tst.rslt2 <- subset(tst.rsltt, grepl("FH", tst.rsltt$W))
    
    cov.tst.rslt11 <- tst.rslt2$Var
    d2$V <- cov.tst.rslt11
    
    d1d2 <- dplyr::full_join(d1,d2, by = c("X1","X2"))
    
    cov.tst.rslt1 <- d1d2$V
    
    cov.tst.1 <- matrix(NA, nrow=length(wt1), ncol=length(wt1))
    
    cov.tst.1[lower.tri(cov.tst.1, diag=FALSE)] <- cov.tst.rslt1
    cov.tst <- t(cov.tst.1)
    cov.tst[lower.tri(cov.tst, diag=FALSE)] <- cov.tst.rslt1
    
    diag(cov.tst) <- var.tst.rslt1
    cov.tst.1 <- matrix(Matrix::nearPD(cov.tst)$mat, length(Z.tst.rslt1),length(Z.tst.rslt1))
    #print(cov.tst.1)
    
    #z.val <- as.vector(ginv(Re(sqrtm(cov.tst)))%*%q.tst.rslt1)
    
    z.max <- max(abs(tst.rslt1$Z))
    cor.tst <- cov2cor(cov.tst.1)
    #print(cor.tst)
    
    #p.value=P(min(Z) < min(z.val))= 1 - P(Z_i >= min(z.val); for all i)
    
    if(one.sided){pval2 <-1 - mvtnorm::pmvnorm(lower = rep(-z.max, length(Z.tst.rslt1)), 
                                               upper= rep(z.max, length(Z.tst.rslt1)), 
                                               corr= cor.tst, algorithm= mvtnorm::GenzBretz())[1]
    
    max.tst <- which(abs(Z.tst.rslt1) == max(abs(Z.tst.rslt1)), arr.ind = TRUE) 
    if(Z.tst.rslt1[max.tst] >= 0){pval <- 1 - pval2/2}  
    if(Z.tst.rslt1[max.tst] < 0){pval <- pval2/2}  
    }
    
    if(!one.sided){pval <- 1 - mvtnorm::pmvnorm(lower = rep(-z.max, length(Z.tst.rslt1)), 
                                                upper= rep(z.max, length(Z.tst.rslt1)), 
                                                corr= cor.tst, 
                                                algorithm= mvtnorm::GenzBretz())[1]}
    
    p.unadjusted <- pnorm(q=tst.rslt1$Z)
    max.index <- which(p.unadjusted == min(p.unadjusted), arr.ind = TRUE)
  }
  
  
  
  #Weighted log-rank test (FH weight): only one weight 
  
  if(!max){
    tst.rslt1 <- subset(tst.rslt, grepl("FH", tst.rslt$W))
    pval <- pnorm(q=tst.rslt1$Z)
    p.unadjusted <- pval
    max.index <- NULL
  }
  
  #Estimation (average HR)
  
  wt.rslt <- data.frame(attr(fit, 'lrw'))
  
  if(max){
    col.wt <- which(grepl("FH", colnames(wt.rslt))==TRUE, arr.ind = T)
    wt.rslt1.1 <- wt.rslt[,c(1,col.wt)]
    wt.rslt1 <- wt.rslt1.1[,max.index]
  }
  
  if(!max){
    col.wt <- which(grepl("FH", colnames(wt.rslt))==TRUE, arr.ind = T)
    wt.rslt1 <- wt.rslt[,col.wt] 
  }
  
  #Performing weighted cox
  
  data.anal.event <- subset(data.anal, status==1)
  data.anal.cens <- subset(data.anal, status==0)
  
  data.anal.event <- data.anal.event[order(data.anal.event$time),]
  
  #Handing ties and matching length of weights
  
  wt.u <- unlist(wt.rslt1)
  
  time.freq <- as.matrix(table(data.anal.event$time))
  
  wt.un <- data.frame(cbind(as.numeric(row.names(time.freq)), time.freq[,1], wt.u))
  wt.all <- wt.un[rep(1:nrow(wt.un), times=wt.un[,2]),]
  data.anal.event$wt <- wt.all$wt.u
  
  #data.anal.event$wt <- unlist(wt.rslt1)
  
  data.anal.cens$wt <- rep(1, nrow(data.anal.cens))
  
  data.anal.w <- rbind(data.anal.event, data.anal.cens)
  data.anal.w$wt[data.anal.w$wt==0] <- 0.000001
  data.anal.w$wt2 <- -log(data.anal.w$wt)
  
  #if(!is.null(strata)){
  
  #  data.anal.w <- cbind(data.anal.w, Strara.mat)
  #  ncol.data.anal.w <- ncol(data.anal.w)
  
  #  data.anal.w <<- data.anal.w
  #  ncol.data.anal.w <<- ncol.data.anal.w
  
  #   FH.est <- coxph(survival::Surv(time, status) ~ arm + strata(data.anal.w[5:ncol.data.anal.w])+ offset(wt2), weight=wt, method= ties.method, data = data.anal.w)
  # }
  
  #if(is.null(strata)){
  
  FH.est <- coxph(survival::Surv(time, status) ~ arm + offset(wt2), weight=wt,  method= ties.method, data = data.anal.w)
  # }
  
  hr <- summary(FH.est)
  hr.est <- hr$conf.int[1]
  hr.low <- hr$conf.int[3]
  hr.up <- hr$conf.int[4]
  
  if(max){
    hr.low.adjusted <- exp(log(hr.est) - (qnorm(1- (alpha/(length(wt) + 1))))*hr$coefficients[3]) 
    hr.up.adjusted <-  exp(log(hr.est) + (qnorm(1- (alpha/(length(wt) + 1))))*hr$coefficients[3]) 
  }
  
  if(max & adjust.methods != "asymp"){out <- list(pval=pval, pval.adjusted = pval.adjusted, p.unadjusted = p.unadjusted, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up, max.index=max.index, hr.low.adjusted= hr.low.adjusted, hr.up.adjusted=hr.up.adjusted,  Z.tst.rslt1= Z.tst.rslt1)}
  if(max & adjust.methods == "asymp"){out <- list(cor=cor.tst,pval=pval, p.unadjusted= p.unadjusted, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up, max.index=max.index, hr.low.adjusted= hr.low.adjusted, hr.up.adjusted=hr.up.adjusted,  Z.tst.rslt1= Z.tst.rslt1, q.tst.rslt1=q.tst.rslt1, max.abs.z=z.max)}
  if(!max){out <- list(pval=pval, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up)}
  
  return(out)
}  
