
summary.cox<-function(data.cox,b.true,rho,gamma){
  bhat<-c(data.cox[,"bhat"])
  hr.true<-exp(b.true)
  avg.hr<-mean(exp(bhat),na.rm=TRUE)
  se.hr<-sqrt(stats::var(exp(bhat),na.rm=TRUE))
  sebhats<-c(data.cox[,"sehat"])
  se.hr.hats<-sqrt((sebhats^2)*((exp(bhat)^2)))
  avg.sehr<-mean(se.hr.hats,na.rm=TRUE)
  hr.lower<-c(data.cox[,"hrL"])
  hr.upper<-c(data.cox[,"hrU"])
  cover.true<-ifelse(c(hr.lower)<=hr.true & c(hr.upper)>=hr.true,1,0)
  cover<-mean(cover.true,na.rm=TRUE)
  
  errl<-ifelse(c(hr.lower)>hr.true,1,0)
  err.low<-mean(errl,na.rm=TRUE)
  
  erru<-ifelse(c(hr.upper)<hr.true,1,0)
  err.up<-mean(erru,na.rm=TRUE)
  
  pow.wald<-mean(c(data.cox[,"R.wald"]),na.rm=TRUE)
  pow.score<-mean(c(data.cox[,"R.score"]),na.rm=TRUE)
  bias<-c(avg.hr-hr.true)
  rel.bias<-bias/hr.true
  ci.length<-mean(hr.upper-hr.lower,na.rm=TRUE)
  cox.results<-as.matrix(cbind(rho,gamma,median(exp(bhat)),
                               avg.hr,se.hr,avg.sehr,cover,err.low,err.up,pow.wald,pow.score))
  colnames(cox.results)<-c("rho","gamma","HR(med)","HR","se(Emp)","se(hat)",
                           "cover","errL","errU",
                           "Pow(Wald)","Pow(Score)")
  return(cox.results)
}




wt.rg.S<-function(S,rho,gamma,tpoints=NULL,t.tau=NULL,w0.tau=0,w1.tau=1){
temp<-S[1:(length(S)-1)]
S.new<-c(1,temp) # S.KM(t-)
term1<-S.new^rho
term2<-(1-S.new)^gamma
wt<-term1*term2
if(!is.null(t.tau)){ 
 if(is.null(tpoints)) stop("tpoints must be provided for specified wt")
   wt<-ifelse(tpoints<=t.tau,w0.tau,w1.tau)  
}
  return(wt)
}


N.rhogamma2<-function(x,error,weight=rep(1,length(error))){
sum(weight*(error<=x))
}

N.rhogamma<-function(x,error,delta,weight=1){
sum(weight*delta*(error<=x))
}


R.rhogamma<-function(x,error,b){
sum(exp(b)*(error>=x))
}

R.rhogamma2<-function(x,error,b,z){
  sum(exp(z*b)*(error>=x))
}

E.rhogamma<-function(x,error,b,z){
  sum(exp(z*b)[which(error==x)])
}



# This version removes the "normalize" option in prior versions
cox.score.rhogamma<-function(beta,time,delta,z,rho,gamma,km.pool,t.tau=NULL,w0.tau=0,w1.tau=1){
is.sorted<-!is.unsorted(time)
if(!is.sorted){
id<-order(time); time<-time[id]; delta<-delta[id]; z<-z[id]
}
at.points<-time

S.pool<-summary(km.pool,c(time))$surv
wt.rg<-wt.rg.S(S=S.pool,rho=rho,gamma=gamma,tpoints=time,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)

# Control group
# Risk and Counting processes
tt<-time[which(z==0)]
dd<-delta[which(z==0)]

risk.z0<-unlist(lapply(as.list(at.points),R.rhogamma,error=tt,b=0.0))
counting<-unlist(lapply(as.list(at.points),N.rhogamma2,error=tt[which(dd==1)]))

N.z0<-counting
dN.z0<-c(wt.rg)*diff(c(0, counting))

tt<-time[which(z==1)]
dd<-delta[which(z==1)]

risk.z1<-unlist(lapply(as.list(at.points),R.rhogamma,error=tt,b=beta))
counting<-unlist(lapply(as.list(at.points),N.rhogamma2,error=tt[which(dd==1)]))
N.z1<-counting
dN.z1<-c(wt.rg)*diff(c(0, counting))

num<-risk.z1*risk.z0
den<-risk.z0+risk.z1
term1<-ifelse(den>0,num/den,0.0)
term2a<-ifelse(risk.z1>0,dN.z1/risk.z1,0.0)
term2b<-ifelse(risk.z0>0,dN.z0/risk.z0,0.0)
score<-sum(term1*(term2b-term2a))
return(score)
}



cox.plik.rhogamma<-function(beta,time,delta,z,rho,gamma,km.pool,t.tau=NULL,w0.tau=0,w1.tau=1){
  is.sorted<-!is.unsorted(time)
  if(!is.sorted){
    id<-order(time); time<-time[id]; delta<-delta[id]; z<-z[id]
  }
  
at.points<-sort(unique(time[which(delta==1)]))

S.pool<-summary(km.pool,c(at.points))$surv
wt.rg<-wt.rg.S(S=S.pool,rho=rho,gamma=gamma,tpoints=at.points,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)

counting<-unlist(lapply(as.list(at.points),N.rhogamma,error=time,delta=delta))
dN<-diff(c(0, counting))
dN.rg<-c(wt.rg)*dN

risk<-unlist(lapply(as.list(at.points),R.rhogamma2,error=time,b=beta,z=z))  
lz<-unlist(lapply(as.list(at.points),E.rhogamma,error=time,b=beta,z=z))  

dN.lik<-dN.rg*c(log(lz)-log(risk))
loglik<-sum(dN.lik)
return(loglik)
}




cox.rhogamma<-function(time,delta,z,rho,gamma,t.tau=NULL,w0.tau=0,w1.tau=1,details=FALSE,est.method="score"){
# x must be a treatment indicator
n<-length(time)
n0<-sum(z==0)
n1<-sum(z==1)
if(n0+n1!=n) stop("x must be a (0/1) treatment group indicator")
is.sorted<-!is.unsorted(time)
if(!is.sorted){
id<-order(time); time<-time[id]; delta<-delta[id]; z<-z[id]
}
# Pooled KM at observed time points for weighting
km.pool<-survfit(Surv(time,delta)~1)


n.ties<-length(time[which(delta==1)])-length(unique(time[which(delta==1)]))

if(est.method=="plik" & n.ties>0){
  cat("Note that PL maximization is not avaiable for weighting in presence of ties:  Setting estimation method to 'score'","\n")
  est.method<-"score"
}

if(est.method=="score"){
  
# Solve score directly
  
get.Cox<-stats::uniroot(f=cox.score.rhogamma,interval=c(-15,15),extendInt="yes",tol=10^-16,
                 time=time,delta=delta,z=z,rho=rho,gamma=gamma,km.pool=km.pool,
                 t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)
if(details){
print(get.Cox)
print(exp(get.Cox$root))
# Check score at 2 points
b2<-jitter(get.Cox$root)
score1<-cox.score.rhogamma(beta=get.Cox$root,time=time,delta=delta,z=z,rho=rho,gamma=gamma,km.pool=km.pool,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)
score2<-cox.score.rhogamma(beta=b2,time=time,delta=delta,z=z,rho=rho,gamma=gamma,km.pool=km.pool,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)
cat("rho,gamma=",c(rho,gamma),"\n")
cat("Score at soln1: beta,score=",c(get.Cox$root,score1),"\n")
cat("Score at b: b,score=",c(b2,score2),"\n")
}
bhat.rhogamma<-get.Cox$root
}

if(est.method=="plik"){

get.Cox<-stats::optimize(f=cox.plik.rhogamma,interval=c(-15,15),maximum=TRUE,
                  time=time,delta=delta,z=z,rho=rho,gamma=gamma,km.pool=km.pool,
                  t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)
if(details){
print(get.Cox)
print(exp(get.Cox$maximum))
}
bhat.rhogamma<-get.Cox$maximum
}

# lik(beta=0)
if(details){
cat("Note: PL calculations are only valid for weighting if NO ties","\n")  
l0<-cox.plik.rhogamma(beta=0,time=time,delta=delta,z=z,rho=rho,gamma=gamma,km.pool=km.pool,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)
l1<-cox.plik.rhogamma(beta=bhat.rhogamma,time=time,delta=delta,z=z,rho=rho,gamma=gamma,km.pool=km.pool,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)
likr<-c(-2*(l0-l1))
cat("rho,gamma=",c(rho,gamma),"\n")
cat("lik0,lik(beta),lik ratio=",c(l0,l1,likr),"\n")
}

# u(beta=0)
u.zero<-cox.score.rhogamma(beta=0,time=time,delta=delta,z=z,rho=rho,gamma=gamma,km.pool=km.pool,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)
u.beta<-cox.score.rhogamma(beta=bhat.rhogamma,time=time,delta=delta,z=z,rho=rho,gamma=gamma,km.pool=km.pool,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)
return(list(bhat=bhat.rhogamma,u.beta=u.beta,u.zero=u.zero))
}





cox.rhogamma.resample<-function(bhat,time,delta,z,G1.draws=NULL,G0.draws=NULL,
                                draws=0,rho=0,gamma=0,details=FALSE,t.tau=NULL,w0.tau=w0.tau,w1.tau=w1.tau,seed.value=8316951){
# The G0.draws and G1.draws will be the same used from max combo test
if(!is.null(G1.draws) & is.null(G0.draws)) stop("G1.draws and G0.draws should be of same dimensions")
if(!is.null(G0.draws) & is.null(G1.draws)) stop("G1.draws and G0.draws should be of same dimensions")
if(!is.null(G0.draws) & !is.null(G1.draws)){ 
if(ncol(G1.draws)!=ncol(G0.draws)) stop("Length of columns should be same (=# of draws)")
}

if(!is.null(G1.draws)) draws<-ncol(G1.draws)

n<-length(time)
n0<-sum(z==0)
n1<-sum(z==1)
if(n0+n1!=n) stop("z must be a (0/1) treatment group indicator")
is.sorted<-!is.unsorted(time)
if(!is.sorted){
id<-order(time); time<-time[id]; delta<-delta[id]; z<-z[id]
}
#at.points<-sort(unique(c(time[delta==1])))
at.points<-time

km.pool<-survfit(Surv(time,delta)~1)
S.pool<-summary(km.pool,c(time))$surv
wt.rg<-wt.rg.S(S=S.pool,rho=rho,gamma=gamma,tpoints=time,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)

tt<-time[which(z==0)]
dd<-delta[which(z==0)]

risk.z0<-unlist(lapply(as.list(at.points),R.rhogamma,error=tt,b=0.0))
counting<-unlist(lapply(as.list(at.points),N.rhogamma,error=tt,delta=dd))
N.z0<-counting
dN.z0<-c(wt.rg)*diff(c(0, counting))

tt<-time[which(z==1)]
dd<-delta[which(z==1)]

risk.z1<-unlist(lapply(as.list(at.points),R.rhogamma,error=tt,b=bhat))
counting<-unlist(lapply(as.list(at.points),N.rhogamma,error=tt,delta=dd))

N.z1<-counting
dN.z1<-c(wt.rg)*diff(c(0, counting))

dN.pooled<-dN.z0+dN.z1

num<-risk.z0*risk.z1
den<-risk.z0+risk.z1
i.bhat<-sum(ifelse(den>0,(num/(den^2))*dN.pooled,0.0))

tt<-time[which(z==0)]
dd<-delta[which(z==0)]

risk.z0<-unlist(lapply(as.list(at.points),R.rhogamma,error=tt,b=0.0))
counting<-unlist(lapply(as.list(at.points),N.rhogamma,error=tt,delta=dd))
N.z0<-counting
dN.z0<- diff(c(0, counting))

tt<-time[which(z==1)]
dd<-delta[which(z==1)]

risk.z1<-unlist(lapply(as.list(at.points),R.rhogamma,error=tt,b=bhat))
counting<-unlist(lapply(as.list(at.points),N.rhogamma,error=tt,delta=dd))
N.z1<-counting
dN.z1<-diff(c(0, counting))

dN.pooled<-dN.z0+dN.z1

Ybar<-risk.z0+risk.z1
Y1<-risk.z1
Y0<-risk.z0
W<-wt.rg
DNbar<-dN.pooled
K<-ifelse(Ybar==0,0,W*(Y1*Y0)/Ybar)


h1<-ifelse(Y1>0,(K^2/Y1),0.0)
h2<-ifelse(Y0>0,(K^2/Y0),0.0)
temp<-Ybar-1
ybar.mod<-ifelse(temp<1,1,temp)
dH1<-ifelse(ybar.mod>0,(DNbar-1)/ybar.mod,0.0)
dH2<-ifelse(Ybar>0,DNbar/Ybar,0.0)
sig2s<-(h1+h2)*(1-dH1)*dH2
sig2.Ub0<-sum(sig2s)

sig.beta<-sqrt(sig2.Ub0/i.bhat^2)
if(details) cat("ASY SE estimates=",c(sig.beta),"\n")

if(draws>0){
set.seed(seed.value)
y0<-time[which(z==0)]
d0<-delta[which(z==0)]

y1<-time[which(z==1)]
d1<-delta[which(z==1)]

n0<-length(y0)
n1<-length(y1)

# (bhat-b0) resamples
bhat.center.star<-rep(NA,draws)
u.bstar<-rep(NA,draws)

if(is.null(G1.draws)){
G0.draws<-matrix(rnorm(draws*n0),ncol=draws)
G1.draws<-matrix(rnorm(draws*n1),ncol=draws)
}

if(nrow(G1.draws)!=n1 | nrow(G0.draws)!=n0 | ncol(G1.draws)!=ncol(G0.draws)) stop("Error in G0.draws and G1.draws")

# U(b0,t) --> perturbed U(bhat,t)

for(dd in 1:draws){
# Z0
g0<-G0.draws[,dd]
# multiply Delta x g0*
counting<-unlist(lapply(as.list(at.points),N.rhogamma,error=y0,delta=d0*g0))
N.z0.star<-counting
dN.z0.star<-diff(c(0, counting))
# Z1
g1<-G1.draws[,dd]
counting<-unlist(lapply(as.list(at.points),N.rhogamma,error=y1,delta=d1*g1))
N.z1.star<-counting
dN.z1.star<-diff(c(0, counting))

num<-W*risk.z1*risk.z0
den<-risk.z0+risk.z1
term1<-ifelse(den>0,num/den,0.0)
term2a<-ifelse(risk.z1>0,dN.z1.star/risk.z1,0.0)
term2b<-ifelse(risk.z0>0,dN.z0.star/risk.z0,0.0)
score.star<-sum(term1*(term2b-term2a))

u.bhat.star<-score.star

u.bstar[dd]<-u.bhat.star
bhat.center.star[dd]<-u.bhat.star/i.bhat

if(details){
if(round(dd/draws,digits=3)==0.10) cat("10% resampling done","\n")
if(round(dd/draws,digits=3)==0.25) cat("25% resampling done","\n")
if(round(dd/draws,digits=3)==0.75) cat("75% resampling done","\n")
if(round(dd/draws,digits=3)==0.90) cat("90% resampling done","\n")
}
}
Z.bstar<-bhat.center.star/sig.beta

var.bhat<-stats::var(bhat.center.star,na.rm=TRUE)
if(details){
#cat("U(bhat) cov. term=",c(sig.z1z0),"\n")
cat("Re-sampling SE estimates=",c(sqrt(var.bhat)),"\n")
#sig.beta<-sqrt(sig2.Ub0/i.bhat^2)
cat("ASY SE estimates=",c(sig.beta),"\n")
}
}
if(draws>0) out<-list(bhat.star=bhat.center.star,var.bhat=var.bhat,u.bstar=u.bstar,i.bhat=i.bhat,sig2.Ub0=sig2.Ub0,se.beta=sig.beta,Z.bstar=Z.bstar)
if(draws==0) out<-list(i.bhat=i.bhat,sig2.Ub0=sig2.Ub0,se.beta=sig.beta)
return(out)
}







fit.cox.rhogamma<-function(time,delta,z,rho,gamma,t.tau=NULL,w0.tau=0,w1.tau=1,draws=0,est.method="score",
                           details=FALSE,seed.value=8316951,G1.draws=NULL,G0.draws=NULL,
                           alphaL=0.025,alphaU=0.025,alphaU.adj=NULL,checking=TRUE){
if(is.null(alphaU.adj)) alphaU.adj<-alphaU
if(!is.null(G1.draws) & is.null(G0.draws)) stop("G1.draws and G0.draws should be of same dimensions")
if(!is.null(G0.draws) & is.null(G1.draws)) stop("G1.draws and G0.draws should be of same dimensions")
if(!is.null(G0.draws) & !is.null(G1.draws)){if (ncol(G1.draws)!=ncol(G0.draws)) stop("Length of columns should be same (=# of draws)")}
  
  
if(!is.null(G1.draws)) draws<-ncol(G1.draws)

fit.rg<-cox.rhogamma(time=time,delta=delta,z=z,rho=rho,gamma=gamma,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau,est.method=est.method)
bhat<-fit.rg$bhat
fit.resample<-cox.rhogamma.resample(bhat=bhat,time=time,delta=delta,z=z,draws=draws,details=details,
                            rho=rho,gamma=gamma,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau,seed.value=seed.value,
                            G1.draws=G1.draws,G0.draws=G0.draws)
bhat<-c(bhat)
se.asy<-fit.resample$se.beta
if(draws>0) se.resample<-sqrt(stats::var(fit.resample$bhat.star))

# Use asymptotic se as default
# For resampling, the bias-corrected is used

sehat<-se.asy

bhatL<-c(bhat+qnorm(alphaL)*sehat)
bhatU<-c(bhat+qnorm(1-alphaU.adj)*sehat)
R.wald<-ifelse(bhatU<0,1,0)

# Bias-correction: BC intervals
# (bhat-b0) distribution
if(draws>0){
b.star<-c(fit.resample$bhat.star+bhat)
j.beta<-mean(b.star<bhat)
z.0<-qnorm(j.beta)
aa<-pnorm(2*z.0+qnorm(alphaL/2))
bb<-pnorm(2*z.0+qnorm(1-alphaU/2))
bhatL.bc<-c(quantile(b.star,c(aa)))
bhatU.bc<-c(quantile(b.star,c(bb)))
R.wald.bc<-ifelse(bhatU.bc<0,1,0)
}
# Score test
u.zero<-fit.rg$u.zero
temp<-cox.rhogamma.resample(bhat=0,time=time,delta=delta,z=z,draws=0,details=FALSE,
                             rho=rho,gamma=gamma,t.tau=t.tau,w0.tau=w0.tau,w1.tau=w1.tau)

sig2.U0<-temp$sig2.Ub0

rm("temp")

z.score<-u.zero/sqrt(sig2.U0)
pvalz<-1-pnorm(z.score)
R.score<-ifelse(pvalz<=alphaU.adj,1,0)


if(rho==0 & gamma==0 & checking){
    # check hr, and ci calculations
  
  checkit<-coxph(Surv(time,delta)~z,ties="breslow")
  ci.check<-c(summary(checkit)$conf.int[c(3,4)])
  hr.check<-exp(checkit$coefficients)
  chisq.check<-checkit$score
  lr0<-survdiff(Surv(time,delta)~z)
if(round(hr.check-exp(bhat),3)!=0){ 
  cat("HR: coxhph vs Mine",c(hr.check,exp(bhat)),"\n")
}
  if(round(ci.check[1]-exp(bhatL),3)!=0){ 
    cat("Lower bound: coxph vs Mine",c(ci.check[1],exp(bhatL)),"\n")
  }
  if(alphaU.adj==0.025 & round(ci.check[2]-exp(bhatU),3)!=0){ 
    cat("Upper Bound: coxph vs Mine",c(ci.check[2],exp(bhatU)),"\n")
  }

  if(round(chisq.check-z.score^2,3)!=0){ 
    cat("Z^2: coxph vs Mine",c(chisq.check,z.score^2),"\n")
    print(lr0)
  }
}  
 
if(details){
U0<-u.zero
z.wald<-bhat/sehat
p.wald<-2*(1-pnorm(abs(z.wald)))
cat("n, num of events=",c(length(time),sum(delta)),"\n")
cat("bhat,se(bhat),z,p(z)=",c(bhat,sehat,z.wald,p.wald),"\n")
cat("hr,hr(lower),hr(upper)=",c(exp(bhat),hr.lower,hr.upper),"\n")
cat("U(0),sig2(U0),Z(lr),p(lr)=",c(U0,sig2.U0,z.score,pvalz),"\n")
cat("Z(lr)^2: Chisq=",c(z.score^2),"\n")
}

if(draws>0){
return(list(rho=rho,gamma=gamma,bhat=bhat,sehat=sehat,se.asy=se.asy,se.resample=se.resample,
bhatL=bhatL,bhatU=bhatU,R.wald=R.wald,
bhatL.bc=bhatL.bc,bhatU.bc=bhatU.bc,R.wald.bc=R.wald.bc,
R.score=R.score,z.score=z.score,pval=pvalz,score=u.zero,sig2=sig2.U0,
hr.hat=exp(bhat),hrL=exp(bhatL),hrU=exp(bhatU),
hrL.bc=exp(bhatL.bc),hrU.bc=exp(bhatU.bc),
fit.resample=fit.resample))
}
if(draws==0){
  return(list(rho=rho,gamma=gamma,bhat=bhat,sehat=sehat,
              bhatL=bhatL,bhatU=bhatU,R.wald=R.wald,
              R.score=R.score,z.score=z.score,pval=pvalz,score=u.zero,sig2=sig2.U0,
              hr.hat=exp(bhat),hrL=exp(bhatL),hrU=exp(bhatU)))
  
}
}





fit.combo.wlr<-function(rgs,time,delta,z,draws,plot.rg=TRUE,print.results=TRUE,outfile=NULL,checking=FALSE,est.method="score",
                        rgsplot=list(c(0,0),c(0,1),c(0,2),c(0,3),c(1,0),c(2,0),c(2,2))){
  if(draws==0) cat("Note that only KM and (rho,gamma)-profile plotting executed: draws>0 necessary for analysis methods")

  if(est.method!="score" & est.method!="plik") stop("Estimation method must be score (solving score for zero root) or plik (max partial-likelihood directly)")
  
    t.start<-proc.time()[1]
  # (rho,gamma) family
  rgs.mat<-matrix(unlist(rgs),ncol=2,byrow=TRUE) # each row one pair (rho,gamma)
  m<-nrow(rgs.mat) # Number of (rho,gamma) combinations
  
  if(plot.rg){
# rgsplot denote the (r,g)'s desired to plot 
    rgs.matplot<-matrix(unlist(rgsplot),ncol=2,byrow=TRUE) # each row one pair (rho,gamma)
    mplot<-nrow(rgs.matplot) # Number of (rho,gamma) combinations
    
    Y1<-time[z==1]; E1<-delta[z==1]
    Y0<-time[z==0]; E0<-delta[z==0]
    
    cls<-c("gray","blue","green","brown","orange","purple","red","orchid")
    
    par(mar = c(5,5,2,5))
    plot(survfit(Surv(Y1,E1)~1),conf.int=FALSE,xlab="time",ylab="Survival",lty=2,col="black",mark.time=TRUE,lwd=4,ylim=c(0,1))
    lines(survfit(Surv(Y0,E0)~1),conf.int=FALSE,lty=2,col="grey",mark.time=TRUE,mark=4,lwd=4)
    
# Note: If rgsplot differs from default, then will need to revise legend (probably best to do outside of this function, then)
# I tried to automate in terms of general (r,g)'s ... but couldn't get legend to work.
        
legend(1,0.6,c("rg=0,1","rg=0,2","rg=0,3","rg=1,0","rg=2,0","rg=2,2"),col=cls[-1],lty=1,lwd=0.5,bty="n")
legend(7,0.9,c("Treatment","Control"),lwd=4,col=c("black","grey"),lty=2,bty="n")
    
    km.pool<-survfit(Surv(Y,Delta)~1)
    atpoints<-sort(unique(Y[Delta==1]))
    S.pool<-summary(km.pool,c(atpoints))$surv
    
    for(rg in 2:mplot){
      
      wt.rg<-wt.rg.S(S=S.pool,rho=rgs.matplot[rg,1],gamma=rgs.matplot[rg,2],tpoints=atpoints)
      
      # scale wt 0-1 range
      wt.rg<-wt.rg/max(wt.rg)
      
      #print(summary(wt.rg))
      
      par(new = TRUE)
      plot(atpoints,wt.rg,axes=FALSE,type="s",lty=1,col=cls[rg],lwd=0.5,xlab=NA,ylab=NA,cex=1.3)
    }
    axis(side = 4)
    mtext(side = 4, line = 3, 'Weights')
  
  } # end plotting of KM along with rho,gamma profiles
  
  # Now, the analyses.
  
  if(draws>0){
  fit.combo<-combo2.logrank(time=time,delta=delta,z=z,draws=draws,
                            rho.gamma=as.list(rgs),
                            include.gehan=FALSE,details=TRUE,
                            checking=checking,draws.status=TRUE)
  
  fit.rgmax<-fit.cox.rhogamma(time=time,delta=delta,z=z,est.method=est.method,
                              rho=fit.combo$rg.max[1],gamma=fit.combo$rg.max[2],                                                                    
                              G1.draws=fit.combo$G1.draws,G0.draws=fit.combo$G0.draws)
  
  # Results to store
  out.names<-c("Z.score","pval(Z)","hr","hrL","hrU","hrL.band","hrU.band","bhat","sehat","se.asy","bhatL","bhatU")
  # Cox model corresponding to max from combo test
  # Resampling
  out.names2<-c("rho","gamma","Zmax","pval(Zmax)","hr","hrL","hrU","hrL.bc","hrU.bc")
  CoxMax.out<-matrix(NA,nrow=1,ncol=length(out.names2))
  colnames(CoxMax.out)<-out.names2
  
  CoxMax.out[1,"rho"]<-c(fit.combo$rg.max)[1]
  CoxMax.out[1,"gamma"]<-c(fit.combo$rg.max)[2]
  
  
  CoxMax.out[1,"Zmax"]<-fit.combo$Z.combo
  CoxMax.out[1,"pval(Zmax)"]<-fit.combo$pval.combo
  
  CoxMax.out[1,"hr"]<-fit.rgmax$hr.hat
  CoxMax.out[1,"hrL"]<-fit.rgmax$hrL
  CoxMax.out[1,"hrU"]<-fit.rgmax$hrU
  CoxMax.out[1,"hrL.bc"]<-fit.rgmax$hrL.bc
  CoxMax.out[1,"hrU.bc"]<-fit.rgmax$hrU.bc
  
  ############################################
  # Initialize for (rho,gamma) combinations
  # Cox.out.1, Cox.out.2, ..., Cox.out.m
  # Cox.out.1 will be standard Cox model
  ############################################
  
  Cox.out.rg<-NULL
  
  for(rg in 1:m){
    temp<-matrix(NA,nrow=1,ncol=length(out.names))
    colnames(temp)<-out.names
    assign(paste("Cox.out",rg,sep="."),temp)
  }
  # Initialize max|Z(r,g)|
  maxZ.beta<-rep(-999,ncol(fit.combo$G1.draws))
  
  for(rg in 1:m){
    fit.rg<-fit.cox.rhogamma(time=time,delta=delta,z=z,rho=rgs.mat[rg,1],gamma=rgs.mat[rg,2],alphaU=0.025,est.method=est.method,
                             G1.draws=fit.combo$G1.draws,G0.draws=fit.combo$G0.draws,checking=checking)
    fit.cox<-fit.rg
    maxzb<-pmax(maxZ.beta,abs(fit.rg$fit.resample$Z.bstar))
    maxZ.beta<-maxzb  # Will update across rg's 
    # temp is assigned as Cox.out.rg
    # Store results for current (r,g)
    # and then rename (update)
    temp<-get(paste("Cox.out",rg,sep="."))
    # Output estimates
    temp[1,"bhat"]<-c(fit.cox$bhat)
    temp[1,"sehat"]<-c(fit.cox$sehat)
    temp[1,"se.asy"]<-c(fit.cox$se.asy)
    temp[1,"bhatL"]<-c(fit.cox$bhatL)
    temp[1,"bhatU"]<-c(fit.cox$bhatU)
    temp[1,"Z.score"]<-c(fit.cox$z.score)  
    temp[1,"pval(Z)"]<-c(fit.cox$pval) 
    temp[1,"hr"]<-c(exp(fit.cox$bhat))
    temp[1,"hrL"]<-c(fit.cox$hrL)
    temp[1,"hrU"]<-c(fit.cox$hrU)
    assign(paste("Cox.out",rg,sep="."),temp)
  }
  
  # Constant for simultaneous band  
  const.rg.band<-quantile(maxZ.beta,c(0.95))
  
  for(rg in 1:m){
    # temp is assigned as Cox.out.rg
    # Store results for current simulation
    # and then rename (update)
    temp<-get(paste("Cox.out",rg,sep="."))
    # Output estimates
    se.asy<-temp[1,"se.asy"]
    bL.band<-temp[1,"bhat"]-const.rg.band*se.asy
    bU.band<-temp[1,"bhat"]+const.rg.band*se.asy
    temp[1,"hrU.band"]<-exp(bU.band)
    temp[1,"hrL.band"]<-exp(bL.band)
    assign(paste("Cox.out",rg,sep="."),temp)
    Cox.out.rg<-rbind(Cox.out.rg,c(rgs.mat[rg,1],rgs.mat[rg,2],temp))
  }
  
  # Focus on WLR Z, and HR estimates
  hr.out<-as.data.frame(Cox.out.rg[,c(1:9)])
  
  colnames(hr.out)<-c("rho","gamma",c(out.names[c(1:7)]))
  
  
  t.end<-proc.time()[1]
  t.min<-(t.end-t.start)/60
  
  
  if(print.results){
    cat("Minutes for resampling draws: draws,mins=",c(draws,t.min),"\n")
    print(CoxMax.out)
    print(hr.out)
  }
  
  if(!is.null(outfile)){ 
    save(CoxMax.out,hr.out,file=outfile)
  }
  
  return(list(CoxMax=CoxMax.out,WLR=hr.out))
} # draws>0
  }






