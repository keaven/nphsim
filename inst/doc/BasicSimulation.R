## ------------------------------------------------------------------------
enrollIntervals <- c(2, 4)

## ------------------------------------------------------------------------
enrollRates <- c(.5, 3, 16)

## ------------------------------------------------------------------------
library(nphsim)
enrollTimes <- rpwexp(n = 300, rate = 7 * enrollRates, intervals = enrollIntervals, cumulative = TRUE)

## ---- fig.cap = "Plot of simulated enrollment."--------------------------
library(ggplot2)
qplot(x = c(0, enrollTimes), y = 0:length(enrollTimes), geom="step", 
      ylab="Number enrolled", xlab="Time") +
      scale_x_continuous(breaks=c(0,6,12,18))

## ------------------------------------------------------------------------
# Failure rates for piecewise exponential time periods
failRates <- c(.3, .6, .5)
# Interval duration(s) before final stable rate
# Note that length is 1 fewer than for failRates
# and should be NULL if there is only 1 failure rate
failIntervals <- c(1,4)

## ------------------------------------------------------------------------
n <- length(enrollTimes)
y <- rpwexp(n = n, rate = failRates, intervals = failIntervals)

## ---- fig.cap = "Plot of sorted simulated failure times without censoring."----
id <- 1:n
# sort patient time-to-event and create a 0 starting point for each patient
dta <- data.frame(N=factor(c(id,id)),Time=c(array(0,n), sort(y, decreasing = TRUE)))
ggplot(dta, aes(x = Time, y = N, grp = N)) + geom_line() +
       xlab("Time-to-event") + ylab("Patients ordered by time-to-event") +
  scale_x_continuous(breaks=(0:4)*6) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())

## ------------------------------------------------------------------------
# Hazard ratio corresponding to control group event rates above
hr <- c(1, .6, .3)
# Sample size of control and experimental arms
ssC <- 5000
ssE <- 5000
# We assume dropouts at a rate of .1 per time unit, increasing to .3 at end
# dropout rates specified for same intervals as failure rates
eta <- c(.1,.1,.3)
# Simulate a single trial instance
# We add an arbitrary enrollment interval for indefinite enrollment duration
trial <- nphsim(nsim = 1, lambdaC = failRates, lambdaE = failRates*hr,
                intervals = failIntervals, ssC = ssC, ssE = ssE,
                gamma = 2 * enrollRates, R = enrollIntervals, fixEnrollTime = FALSE,
                eta = eta, etaE = eta)
# show a few lines
head(trial$simd, n=5)

## ------------------------------------------------------------------------
table(trial$simd$treatment, trial$simd$cnsr)

## ------------------------------------------------------------------------
summary(trial$simd$enterT+trial$simd$survival)

## ------------------------------------------------------------------------
library(survival)
plot(with(trial$simd, survfit(Surv(survival, 1-cnsr) ~ treatment)))

## ------------------------------------------------------------------------
testOut <- simtest(x = trial, anaD=900, method='LR')
testOut$result

## ----message=FALSE, warning=FALSE----------------------------------------
testOut <- simtest(x = trial, anaD=900, method=wlr.Stat, fparam=list(rho=c(0,1), gamma=c(1,1)))
testOut$result

