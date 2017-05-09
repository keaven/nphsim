nSurvNPH<-function ( lambdaC=log(2)/6
                    ,hr = 0.6
                    ,hr0 = 1
                    ,etaC = 0
                    ,etaE = NULL
                    ,NONCMPL=0
                    ,DROPIN= 0
                    ,gamma = 1
                    ,R = 12
                    ,S = NULL
                    ,T = NULL
                    ,minfup = NULL
                    ,ratio = 1
                    ,alpha = 0.025
                    ,beta = 0.1
                    ,sided = 1, tol = .Machine$double.eps^0.25
                    ,SBDV=20
                    ,SIMULT=FALSE
                    ,TRT_LAG=0 #lag treatment#
                    ,nstrata=1)

{
  zalpha <- -qnorm(alpha/sided)
  zbeta<--qnorm(beta)

  Qe<-ratio/(1+ratio) #proprotion of Experimental group
  Qc<-1-Qe #proportion of Control group

  N_INTRVL <- SBDV
  if (TRT_LAG>0)
  { a<-SBDV*TRT_LAG
    if (abs(a-round(a))>1e-6) stop("SBDV*TRT_LAG should be integer")
    NACTV    <- round(a)
    NSTATES  <- (2 * NACTV) + 2;
  }

  NN<-round(T*N_INTRVL) #total number of time points in Markov process
  #adjust gamma and recrutiment period given R and T;
  if (!is.matrix(gamma)) gamma <- matrix(gamma)
  cumR_org<-cumsum(R)
  cumR_adj<-unique(c(0,cumR_org[cumR_org<T],min(T,sum(R))))
  R<-diff(cumR_adj)
  nper<-length(R) #no. of recruitment period
  gamma<-matrix(gamma[1:nper,],ncol=nstrata)
  mat0<-matrix(0,NN,nstrata)
  gamma_NN<-mat0
  temp<-apply(gamma,2,FUN=rep,times=R*N_INTRVL)/N_INTRVL
  gamma_NN[1:nrow(temp),]<-temp
  RCRT_CUM<-matrix(apply(gamma_NN,2,cumsum),ncol=nstrata)
  AD_CENS<-mat0
  AD_CENS[NN:1,]<-gamma_NN/RCRT_CUM; #administrative censoring;

  if (is.null(S)) {nper<-1;S<-T}
  else nper<-length(S)+1
  if (!is.matrix(lambdaC)) lambdaC <- matrix(lambdaC, nrow=nper, ncol=nstrata) #row represent time, col represent strata
  if (!is.matrix(NONCMPL)) NONCMPL <- matrix(NONCMPL, nrow=nper, ncol=nstrata)
  if (!is.matrix(DROPIN)) DROPIN <- matrix(DROPIN, nrow=nper, ncol=nstrata)
  if (!is.matrix(hr)) hr <- matrix(hr,nrow=nper,ncol=nstrata)
  if (is.null(etaE)) etaE<-etaC
  if (!is.matrix(etaC)) etaC <- matrix(etaC,nrow=nper,ncol=nstrata)
  if (!is.matrix(etaE)) etaE <- matrix(etaE,nrow=nper,ncol=nstrata)

  #adjust event rate, drop out rate given S and T;
  cumS_org<-cumsum(S)
  cumS_adj<-unique(c(0,cumS_org[cumS_org<T],T))
  S<-diff(cumS_adj)
  lambdaE<-hr*lambdaC

  SN<-S*N_INTRVL
  nper<-length(S)

  NONCMPL_NP<-1-(1-NONCMPL[1:nper,,drop=FALSE])^(1/(N_INTRVL))
  DROPIN_NP<-1-(1-DROPIN[1:nper,,drop=FALSE])^(1/(N_INTRVL))
  etaC_NP<-1-(1-etaC[1:nper,,drop=FALSE])^(1/(N_INTRVL))
  etaE_NP<-1-(1-etaE[1:nper,,drop=FALSE])^(1/(N_INTRVL))
  pC_NP<-1-exp(-lambdaC[1:nper,,drop=FALSE])^(1/N_INTRVL)
  pE_NP<-1-exp(-lambdaE[1:nper,,drop=FALSE])^(1/N_INTRVL)


  NONCMPL_NN<-apply(NONCMPL_NP,2,FUN=rep,times=SN)
  DROPIN_NN<- apply(DROPIN_NP,2,FUN=rep,times=SN)
  etaC_NN<-apply(etaC_NP,2,FUN=rep,times=SN)#drop out rate for control;
  etaE_NN<-apply(etaE_NP,2,FUN=rep,times=SN)#drop out rate for experiment;
  pC_NN<-  apply(pC_NP,2,FUN=rep,times=SN)  #event rate for control;
  pE_NN<-  apply(pE_NP,2,FUN=rep,times=SN)  #event rate for experiment;

#MARKOV;
#INITAILIZE MATRICES;
DSTR<-vector("list", nstrata)


if (!TRT_LAG)
{ P_C<-rep(0,nstrata)
  P_E<-rep(0,nstrata)
  P_ALL<-rep(0,nstrata)
  eD_LR<-rep(0,nstrata)
  for (strata in 1:nstrata)
  {
    #initial state
  DISTR_E<-matrix(c(0,0,1,0),ncol=1)
  DISTR_C<-matrix(c(0,0,0,1),ncol=1)
  TRANS<-matrix(0,4,4)
  DSTR_E<-matrix(0,4,NN)
  DSTR_C<-matrix(0,4,NN)
  rownames(DSTR_C)<-c('LOSSES', 'EVENTS', 'DROPIN', 'ACTV C')
  rownames(DSTR_E)<-c('LOSSES', 'EVENTS', 'ACTV_E', 'NONCMPL')

   #transition matrix
    for ( I in 1:NN)
    { TRANS<-diag(4);
      TRANS[,3]<-c(etaE_NN[I,strata],pE_NN[I,strata],0,NONCMPL_NN[I,strata]);
      TRANS[,4]<-c(etaC_NN[I,strata],pC_NN[I,strata],DROPIN_NN[I,strata],0);
      diag(TRANS)<-0
      diag(TRANS)<-1-colSums(TRANS);
      DISTR_E<-TRANS%*%DISTR_E;
      DISTR_C<-TRANS%*%DISTR_C;

      if (!SIMULT) #start addressing stagged entry
        {
        TEMP_E<-DISTR_E[c(3,4)]*(1-AD_CENS[I,strata]);
        DISTR_E[1]<-DISTR_E[1]+sum(DISTR_E[c(3, 4)]-TEMP_E);
        DISTR_E[c(3,4)]<-TEMP_E;
        TEMP_C<-DISTR_C[c(3,4)]*(1-AD_CENS[I,strata]);
        DISTR_C[1]<-DISTR_C[1]+sum(DISTR_C[c(3, 4)]-TEMP_C);
        DISTR_C[c(3,4)]<-TEMP_C;
        }#end of addressing stagged entry


      DSTR_E[,I]<-DISTR_E;
      DSTR_C[,I]<-DISTR_C;
      }

  EVENT_C<-diff(c(0,DSTR_C[2,]))
  EVENT_E<-diff(c(0,DSTR_E[2,]))
  EVENT_ALL<-Qc*EVENT_C+Qe*EVENT_E
  LOSS_C<-diff(c(0,DSTR_C[1,]))
  LOSS_E<-diff(c(0,DSTR_C[1,]))
  ATRISK_C<-1-DSTR_C[1,]-DSTR_C[2,]+EVENT_C+LOSS_C
  ATRISK_E<-1-DSTR_E[1,]-DSTR_E[2,]+EVENT_E+LOSS_E
  PHI<-Qc/Qe*ATRISK_C/ATRISK_E #population ratio
  THETA<-log(1-EVENT_C/ATRISK_C)/log(1-EVENT_E/ATRISK_E) #risk ratio
  RHO<-EVENT_ALL/sum(EVENT_ALL)  # di/d
  GAMMA<-PHI*THETA/(1+PHI*THETA)-PHI/(1+PHI);
  ETA<-PHI/(1+PHI)^2;
  P_E[strata]<-DSTR_E[2,ncol(DSTR_E)]
  P_C[strata]<-DSTR_C[2,ncol(DSTR_C)]
  P_ALL[strata]<-Qc*P_C[strata]+Qe*P_E[strata]
  eD_LR[strata]<-sum(RHO*GAMMA)/sqrt(sum(RHO*ETA))
  DSTR[[strata]]<-list(DSTR_E=DSTR_E, DSTR_C=DSTR_C)

  }

} #end of no treatment lag

else if (TRT_LAG>0)
  {
  lambdaC<-lambdaC[1,1,drop=FALSE]
  lambdaE<-lambdaE[1,1,drop=FALSE]
  P_C<-rep(0,nstrata)
  P_E<-rep(0,nstrata)
  P_ALL<-rep(0,nstrata)
  eD_LR<-rep(0,nstrata)
  for (strata in 1:nstrata)
  {
    A<-matrix(0,NACTV,NACTV)
    B<-matrix(0,NACTV,NACTV)
    C<-matrix(0,NACTV,NACTV)
    D<-matrix(0,NACTV,NACTV)

    #initial state
    DISTR_E<-matrix(0,nrow=NSTATES,ncol=1)
    DISTR_C<-matrix(0,nrow=NSTATES,ncol=1)
    DISTR_E[3]<-1
    DISTR_C[NACTV+3]<-1
    TRANS<-matrix(0,NSTATES,NSTATES)
    DSTR_E<-matrix(0,4,NN)
    DSTR_C<-matrix(0,4,NN)
    rownames(DSTR_C)<-c('LOSSES', 'EVENTS', 'DROPIN', 'ACTV_C')
    rownames(DSTR_E)<-c('LOSSES', 'EVENTS', 'ACTV_E', 'NONCMPL')


    #transition matrix
    for ( I in 1:NN)
    {TRANS<-matrix(0,NSTATES,NSTATES);
    A<-matrix(0,NACTV,NACTV)
    B<-matrix(0,NACTV,NACTV)
    C<-matrix(0,NACTV,NACTV)
    D<-matrix(0,NACTV,NACTV)
    AA<-1-exp(-lambdaC-(0:NACTV)*(lambdaE-lambdaC)/NACTV)^(1/N_INTRVL)
    TRANS[1,]<-c(1,0,rep(etaE_NN[I,strata],NACTV),rep(etaC_NN[I,strata],NACTV))
    TRANS[2,]<-c(0,1,AA[2:(NACTV+1)],AA[1:(NACTV)])

    C<-diag(rep(DROPIN_NN[I,strata],NACTV))
    B<-diag(rep(NONCMPL_NN[I,strata],NACTV))

    if (NACTV>2) {
    diag(A[2:NACTV,1:(NACTV-1)])<-1-colSums(TRANS[,3:(NACTV+1),drop=FALSE])-NONCMPL_NN[I,strata]
    diag(D[1:(NACTV-1),2:(NACTV)])<-1-colSums(TRANS[,(NACTV+4):(2*NACTV+2),drop=FALSE])-DROPIN_NN[I,strata]
    }
    else
    {A[2,1]<-1-sum(TRANS[,3])-NONCMPL_NN[I,strata]
     D[1,2]<-1-sum(TRANS[,6])-DROPIN_NN[I,strata]
     }
    A[NACTV,NACTV]<-1-colSums(TRANS[,(NACTV+2),drop=FALSE])-NONCMPL_NN[I,strata]


    D[1,1]<-1-colSums(TRANS[,(NACTV+3),drop=FALSE])-DROPIN_NN[I,strata]
    TRANS[(3:NSTATES),(3:NSTATES)]<-cbind(rbind(A,B),rbind(C,D))
    DISTR_E<-TRANS%*%DISTR_E;
    DISTR_C<-TRANS%*%DISTR_C;

    if (!SIMULT) #start addressing stagged entry
    {
      TEMP_E<-DISTR_E[c(3:NSTATES)]*(1-AD_CENS[I,strata]);
      DISTR_E[1]<-DISTR_E[1]+sum(DISTR_E[c(3:NSTATES)]-TEMP_E);
      DISTR_E[3:NSTATES]<-TEMP_E;
      TEMP_C<-DISTR_C[3:NSTATES]*(1-AD_CENS[I,strata]);
      DISTR_C[1]<-DISTR_C[1]+sum(DISTR_C[3:NSTATES]-TEMP_C);
      DISTR_C[3:NSTATES]<-TEMP_C;
    }#end of addressing stagged entry


    DSTR_E[,I]<-c(DISTR_E[1:2],sum(DISTR_E[3:(NACTV+2)]),sum(DISTR_E[(NACTV+3):NSTATES]))
    DSTR_C[,I]<-c(DISTR_C[1:2],sum(DISTR_C[3:(NACTV+2)]),sum(DISTR_C[(NACTV+3):NSTATES]))
    }

    EVENT_C<-diff(c(0,DSTR_C[2,]))
    EVENT_E<-diff(c(0,DSTR_E[2,]))
    EVENT_ALL<-Qc*EVENT_C+Qe*EVENT_E
    LOSS_C<-diff(c(0,DSTR_C[1,]))
    LOSS_E<-diff(c(0,DSTR_C[1,]))
    ATRISK_C<-1-DSTR_C[1,]-DSTR_C[2,]+EVENT_C+LOSS_C
    ATRISK_E<-1-DSTR_E[1,]-DSTR_E[2,]+EVENT_E+LOSS_E
    PHI<-Qc/Qe*ATRISK_C/ATRISK_E #population ratio
    THETA<-log(1-EVENT_C/ATRISK_C)/log(1-EVENT_E/ATRISK_E) #risk ratio
    RHO<-EVENT_ALL/sum(EVENT_ALL)  # di/d
    GAMMA<-PHI*THETA/(1+PHI*THETA)-PHI/(1+PHI);
    ETA<-PHI/(1+PHI)^2;
    P_E[strata]<-DSTR_E[2,ncol(DSTR_E)]
    P_C[strata]<-DSTR_C[2,ncol(DSTR_C)]
    P_ALL[strata]<-Qc*P_C[strata]+Qe*P_E[strata]
    eD_LR[strata]<-sum(RHO*GAMMA)/sqrt(sum(RHO*ETA))
    DSTR[[strata]]<-list(DSTR_E=DSTR_E, DSTR_C=DSTR_C)

  }




  }  #end of treatment lag
q_strata<-colSums(gamma_NN)/sum(gamma_NN)
w_strata<-q_strata*P_E*P_C/P_ALL
w_strata<-w_strata/sum(w_strata)
dj<-q_strata*P_ALL/sum(q_strata*P_ALL)#proportion death for each stratum
C1<-sum(w_strata*eD_LR*dj^0.5)
d<-(zalpha+zbeta)^2/C1^2*sum(w_strata^2)
n<-d/sum(q_strata*P_ALL)

eNC<-n*q_strata*Qc
eNE<-n*q_strata*Qe
eN<-n*q_strata
eDC<-eNC*P_C
eDE<-eNE*P_E
eD<-eN*P_ALL

output<-list(DSTR=DSTR,alpha=alpha, beta=beta
             , eDC=eDC,eDE=eDE,d=eD
             , eNC=eNC,eNE=eNE,n=eN
             , pE=P_E, pC=P_C, p=P_ALL
             , d=round(d),n=round(n))
return(output)

}


#nSurvNPH();
#no staggered entry, no lag time;

#no staggered entry, lag time;
#stagger entry, no lag time
#stagged entry; lag time;
#strata;
#undebug("nSurvNPH")
#source("nSurvNPHDebug.R")


