cox_pw <- function(time     = NULL,
                   status   = NULL,
                   arm      = NULL,
                   cutpt    = NULL
)
  
{   
  data <- data.frame(time,status, arm)
  colnames(data) <- c("time", "status", "arm")
  
  #data$id <- 1:length(time)
  
  dd <-survSplit(data,cut= cutpt,end="time",start="start",event="status", id="id", episode = "tint")
  
  coxrslt <- coxph(Surv(start, time, status) ~ arm:strata(tint), data=dd)
  
  hr <- summary(coxrslt)
  est <- data.frame(cbind(hr$conf.int[,1], hr$conf.int[,2], hr$conf.int[,3], hr$conf.int[,4]))
  colnames(est) <- c("HR", "SE", "95% Lower", "95% Upper")
  
  nn <- length(unique(dd$tint))
  row.names(est) <- paste("treatment (", "interval", 1:nn, ")", sep='')
  
  
  return(list(est=est, data=data))
  
}
