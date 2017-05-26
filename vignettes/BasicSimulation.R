## ------------------------------------------------------------------------
enrollIntervals <- c(2, 4, 12)

## ------------------------------------------------------------------------
enrollRates <- c(1, 2, 4)

## ------------------------------------------------------------------------
library(nphsim)
enrollTimes <- rpwexp(n = 300, rate = 10 * enrollRates, intervals = enrollIntervals, cumulative = TRUE)

## ---- fig.cap = "Plot of simulated enrollment."--------------------------
library(ggplot2)
qplot(x = c(0, enrollTimes), y = 0:length(enrollTimes), geom="step", 
      ylab="Number enrolled", xlab="Time") +
      scale_x_continuous(breaks=c(0,2,6,10))

## ------------------------------------------------------------------------
# Failure rates for piecewise exponential time periods
failRates <- c(.8, .4, .2, .1)
# Interval durations before final stable rate
# Note that length is 1 fewer than for failRates
failIntervals <- c(0, 3, 6)

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
# and corresponding failure rates for experimental group
hr <- c(1, .7, .5, .3)
lambdaE <- hr * failRates
# Sample size of control and experimental arms
ssC <- n
ssE <- n
# We assume dropouts at a rate of .05, increasing to .1 at end
# dropout rates specified for same intervals as failure rates
eta <- c(.1,.1,.1,.3)
# Simulate a single trial instance
# We add an arbitrary enrollment interval for indefinite enrollment duration
trial <- nphsim(nsim = 1, lambdaC = failRates, lambdaE = lambdaE,
                intervals = failIntervals, ssC = ssC, ssE = ssE,
                gamma = 2 * enrollRates, R = c(enrollIntervals),
                eta = eta, etaE = eta)
# show a few lines
head(trial$simd, n=5)

## ------------------------------------------------------------------------
table(trial$simd$treatment, trial$simd$cnsr)

## ------------------------------------------------------------------------
summary(trial$simd$enterT+trial$simd$survival)

## ------------------------------------------------------------------------
library(survival)
plot(with(trial$simd, survfit(Surv(survival, 1- cnsr) ~ treatment)))

## ------------------------------------------------------------------------
testOut <- simtest(x = trial, anaD=300, method='LR')
testOut$result

