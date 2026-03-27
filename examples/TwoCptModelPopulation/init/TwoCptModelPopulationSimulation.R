## Template to simulate PKPD data for population model
## Inter-individual variability on all PK parammeters

rm(list = ls())
gc()

modelName <- "TwoCptModelPopulation"

library(rstan)
library(mrgsolve)

###############################################################################################
## mrgsolve

nSub <- 10; # number of subjects
nIIV <- 5; # number of parameters with Inter-Individual variation

code <- '
$PARAM CL = 5, Q = 8, V2 = 20, V3 = 70, KA = 1.2

$SET delta = 0.1  // simulation grid

$CMT GUT CENT PERI

$GLOBAL
#define CP (CENT/V2)

$PKMODEL ncmt = 2, depot = TRUE

$MAIN
double CLi = exp(log(CL) + ETA(1));
double Qi = exp(log(Q) + ETA(2));
double V2i = exp(log(V2) + ETA(3));
double V3i = exp(log(V3) + ETA(4));
double KAi = exp(log(KA) + ETA(5));

$OMEGA  name="IIV"
0.0025 0.0025 0.0025 0.0025 0.0025

$SIGMA 0.01

$TABLE
capture DV = CP * exp(EPS(1));

$CAPTURE DV CP
'

mod <- mread("accum", tempdir(), code) %>% Req(DV) %>% update(end=480, delta=0.1)
e1 <- expand.ev(amt = rep(80 * 1000, nSub), ii = 12, addl = 14) # Create dosing events
out <- mod %>% data_set(e1) %>% carry.out(dose) %>% Req(DV) %>% mrgsim(end = 250)
plot(out, DV~time|factor(ID),scales="same")

## Observation times
time <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
time <- c(time, time + 12, seq(24, 156, by = 12), c(time, 12, 18, 24) + 168)
time <- time[time != 0]

# save data in data frame 
SimData <- 
  mod %>%
  data_set(e1) %>%
  carry.out(cmt, ii, addl, rate, amt, evid, ss) %>%
  mrgsim(Req = "DV", end = -1, add = time, recsort = 3) %>%
  as.data.frame

head(SimData)

SimData$cmt[SimData$cmt == 0] <- 2 ## adjust cmt (adopt NONMEM convention)
## remove observation with 0 drug concentration
SimData <- SimData[!((SimData$evid == 0)&(SimData$DV == 0)),]

################################################################################################
## Format data for Stan 
nt <- nrow(SimData)
iObs <- with(SimData, (1:nrow(SimData))[evid == 0])
nObs <- length(iObs)

## Subject specific data
xsub <- subset(SimData, !duplicated(ID))
nSubjects <- length(xsub$ID)

## Row indices for start and end of each individual's data
start <- (1:nt)[!duplicated(SimData$ID)]
end <- c(start[-1] - 1, nt)

## create Stan data set
data <- with(SimData,
             list(nt = nt,
                  nSubjects = nSubjects,
                  start = start,
                  end = end,
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
                  ss = ss,
                  nIIV = nIIV)) ## number of parameters with IIV

## create initial estimates
init <- function()
  list(CLHat = exp(rnorm(1, log(10), 0.2)),
       QHat = exp(rnorm(1, log(20), 0.2)),
       V1Hat = exp(rnorm(1, log(70), 0.2)),
       V2Hat = exp(rnorm(1, log(70), 0.2)),
       kaHat = exp(rnorm(1, log(1), 0.2)),
       sigma = runif(1, 0.5, 2),
       L = diag(nIIV),
       etaStd = matrix(rep(0, nIIV * data$nSubjects), nrow = nIIV),
       omega = runif(nIIV, 0.5, 2),
       logtheta = matrix(rep(log(c(exp(rnorm(1, log(10), 0.2)),
                                   exp(rnorm(1, log(20), 0.2)),
                                   exp(rnorm(1, log(70), 0.2)),
                                   exp(rnorm(1, log(70), 0.2)),
                                   exp(rnorm(1, log(1), 0.2)))),
                             ea = data$nSubjects),
                         nrow = data$nSubjects))

with(data, stan_rdump(ls(data), file = paste0(modelName,".data.R")))
inits <- init()
with(inits, stan_rdump(ls(inits), file = paste0(modelName,".init.R")))
