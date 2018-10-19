## Template to simulate PKPD data
rm(list = ls())
gc()

modelName <- "TwoCptModel"

library(rstan)
library(mrgsolve)  ## tested with version 0.7.6.9029

set.seed(11091989) ## not required but assures repeatable results

###############################################################################################
## Simulate data using mrgsolve

code <- '
$PARAM CL = 5, Q = 8, V2 = 20, V3 = 70, KA = 1.2

$CMT GUT CENT PERI

$GLOBAL
#define CP (CENT/V2)
 
$PKMODEL ncmt = 2, depot = TRUE

$SIGMA 0.01  // variance

$TABLE
capture DV = CP * exp(EPS(1));

$CAPTURE CP
'

mod <- mcode("accum", code) %>% Req(DV) %>% update(end=480,delta=0.1)

e1 <- ev(amt = 80 * 1000, ii = 12, addl = 14) # Create dosing events
mod %>% ev(e1) %>% mrgsim(end = 250) %>% plot # plot data

## Observation times
time <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
time <- c(time, time + 12, seq(24, 156, by = 12), c(time, 12, 18, 24) + 168)
time <- time[time != 0]

# save data in data frame 
SimData <- 
  mod %>% 
  ev(e1) %>% 
  carry_out(cmt, ii, addl, rate, amt, evid, ss) %>%
  mrgsim(Req = "DV", end = -1, add = time, recsort = 3) %>%
  as.data.frame

head(SimData)

SimData$cmt[SimData$cmt == 0] <- 2 ## adjust cmt (adopt NONMEM convention)
SimData <- SimData[!((SimData$evid == 0)&(SimData$DV == 0)),] ## remove observation with 0 drug concentration

################################################################################################
# Format Data for Stan using RStan
nt <- nrow(SimData)
iObs <- with(SimData, (1:nrow(SimData))[evid == 0])
nObs <- length(iObs)

## create Stan data set
data <- with(SimData,
             list(nt = nt,
                  nObs = nObs,
                  iObs = iObs,
                  time = time,
                  cObs = DV[iObs],
                  amt =  amt,
                  rate = rate,
                  cmt = cmt,
                  evid = evid,
                  ii = ii,
                  addl = addl,
                  ss = ss))

## create initial estimates
## Note: you really want to create different initial estimates for each
## chain. This is done TwoCptCmdStan.R.
init <- function() 
  list(CL = exp(rnorm(1, log(10), 0.2)),
       Q = exp(rnorm(1, log(20), 0.2)),
       V1 = exp(rnorm(1, log(70), 0.2)),
       V2 = exp(rnorm(1, log(70), 0.2)),
       ka = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2))

with(data, stan_rdump(ls(data), file = paste0(modelName, ".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName, ".init.R")))
