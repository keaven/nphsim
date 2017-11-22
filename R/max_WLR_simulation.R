

rm.combo.WLRmax <- function(time        = NULL,
                            status      = NULL,
                            arm         = NULL,
                            strata      = NULL, #Not working for p-value at this moment due to SurMisc package does not support StrataTen
                            wt          = NULL,
                            ties.method = c("exact", "breslow", "efron")[2],
                            adjust.methods = c("holm", "hochberg", "hommel", "bonferroni")[3],
                            max         = TRUE,
                            alpha       = 0.05
)
{
  
  if(!is.null(strata)){
    
  n.obs <- sapply(strata,length)
  seq.max <- seq_len(max(n.obs))
  Strara.mat <- sapply(strata, "[", i = seq.max)

  #data.anal <- data.frame(cbind(time,status,arm, Strara.mat))
  #ncol.data.anal <- ncol(data.anal)
  #data.anal <<- data.anal
  #ncol.data.anal <<- ncol.data.anal
  #fit<- ten(survfit(Surv(time, status) ~ arm + strata(data.anal[,4:ncol.data.anal]), data = data.anal))  
                      }
  
  #if(is.null(strata)){
    data.anal <- data.frame(cbind(time,status,arm))
    fit<- ten(survfit(Surv(time, status) ~ arm, data = data.anal))  
   #                  }

#Testing
  
  comp(fit, p= sapply(wt, function(x){x[1]}), q= sapply(wt, function(x){x[2]}))
  
  tst.rslt <- attr(fit, 'lrt')

#Combination test (Using Holms-Bonferroni for muttiplicity adjustment)
  
  if(max){
    
    tst.rslt1 <- rbind(tst.rslt[1,],subset(tst.rslt, grepl("FH", tst.rslt$W)))
    max.tst.rslt1 <- max(tst.rslt1$Z)
    pval.unadjusted <- pnorm(q=tst.rslt1$Z)
    #print(pval.unadjusted)
    #pval.adjusted <- p.adjust(pval.unadjusted, method = "holm")
    pval.adjusted <- p.adjust(pval.unadjusted, method = adjust.methods)
    pval <- min(pval.adjusted)
    #max.index <- which(pval == min(pval.adjusted), arr.ind = TRUE)
    max.index <- which(pval.unadjusted == min(pval.unadjusted), arr.ind = TRUE)
    #print(max.index)
  
  }
  
  
#Weighted log-rank test (FH weight)
  
  if(!max){
    
    tst.rslt1 <- subset(tst.rslt, grepl("FH", tst.rslt$W))
    pval <- pnorm(q=tst.rslt1$Z)
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
  
 
  if(!is.null(strata)){
  
  data.anal.w <- cbind(data.anal.w, Strara.mat)
  ncol.data.anal.w <- ncol(data.anal.w)
  
  data.anal.w <<- data.anal.w
  ncol.data.anal.w <<- ncol.data.anal.w

  FH.est <- coxph(Surv(time, status) ~ arm + strata(data.anal.w[5:ncol.data.anal.w])+ offset(wt2), weight=wt, method= ties.method, data = data.anal.w)
  
  }
  
  if(is.null(strata)){
    
    FH.est <- coxph(Surv(time, status) ~ arm + offset(wt2), weight=wt,  method= ties.method, data = data.anal.w)
    
                     }
  
  
  hr <- summary(FH.est)
  hr.est <- hr$conf.int[1]
  hr.low <- hr$conf.int[3]
  hr.up <- hr$conf.int[4]
  
  
  if(max){
    hr.low.adjusted <- exp(log(hr.est) - (qnorm(1- (alpha/(length(wt) + 1))))*hr$coefficients[3]) 
    hr.up.adjusted <-  exp(log(hr.est) + (qnorm(1- (alpha/(length(wt) + 1))))*hr$coefficients[3]) 
    
  }
  
  if(max){out <- list(pval=pval, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up, max.index=max.index, hr.low.adjusted= hr.low.adjusted, hr.up.adjusted=hr.up.adjusted)}
  if(!max){out <- list(pval=pval, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up)}
  
  return(out)
    
}  

# 
# nsim <- 10
# medC <- c(6,6)
# hr <- c(1,1)
# intervals <- c(10) # intervals for the piecewise failure rate is set up per HR defined above
# gamma <- c(2.5, 5) ## a ramp-up enrollment
# R     <- c(2, 12) ## enrollment period: total of 12 months
# eta <- -log(0.90) ## 1% monthly dropout rate
# 
# sim10 <- nphsim(nsim=nsim, lambdaC=log(2)/medC, lambdaE=(log(2)/medC)*hr,
#                 ssC=300, ssE=300, intervals=intervals, gamma=gamma, R=R, eta=eta)



# ###Testing the macro
# 
# tt <- sim10$simd[sim==3,]$survival
# trt <- sim10$simd[sim==3,]$treatment
# cens <- 1-sim10$simd[sim==3,]$cnsr
# strata <- list(st1=rbinom(600,1,0.5), st2=rbinom(600,1,0.4))

# #FH(0,1)
# rm.combo.WLRmax(time= tt, status = cens, arm= trt, wt= list(a1=c(0,1)), max= FALSE, strata=strata)
# 
#Combo of FH(1,0), FH(1,1), FH(0,1)

#rm.combo.WLRmax(time= tt, status = cens, arm= trt, wt= list(a1=c(0,1), a2=c(1,1), a3=c(1,0)), max= TRUE)

test.combo <- function(survival, cnsr, treatment, stratum=NULL, fparam) {
  test.t <- rm.combo.WLRmax(time=survival, status=1-cnsr, arm=(treatment=='experimental'), wt=fparam, max=T)
  pval <- test.t$pval
  hr.est <- test.t$hr.est
  hr.low <- test.t$hr.low 
  hr.up <- test.t$hr.up 
  max.index <- test.t$max.index
  
  return(list(pval=pval, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up, max.index=max.index))
}  

# ms.test <- simtest(x=sim10, anaD=450, method=test.combo, stratum=NULL, fparam=list(a1=c(0,1), a2=c(1,1), a3=c(0,3)))
# ms.test$result


test.WLR <- function(survival, cnsr, treatment, stratum=NULL, fparam) {
  test.t <- rm.combo.WLRmax(time=survival, status=1-cnsr, arm=(treatment=='experimental'), wt=fparam, max=F)
  pval <- test.t$pval
  hr.est <- test.t$hr.est
  hr.low <- test.t$hr.low 
  hr.up <- test.t$hr.up 
  #max.index <- test.t$max.index
  
  return(list(pval=pval, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up))
} 


