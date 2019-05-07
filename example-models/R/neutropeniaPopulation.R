rm(list = ls())
gc()

setwd("~/Desktop/Code/torsten/example-models/PKPD/torsten/R")
.libPaths("~/svn-StanPmetrics/script/lib")
modelName <- "neutropeniaPopulation"

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path(scriptDir, "deliv", "figure", modelName)
tabDir <- file.path(scriptDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, modelName)
modelDir <- file.path(projectDir, modelName)
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "tools")

stanDir <- file.path(scriptDir, "cmdstan")
tempDir <- file.path(scriptDir, "temp")
dir.create(tempDir)

library(rstan)
library(metrumrg)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(loo)
library(parallel)
library(mrgsolve)
source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "functions.R"))

rstan_options(auto_write = TRUE)

set.seed(11191951) ## not required but assures repeatable results

################################################################################################
## Define function to run model with fixed parameters.

runModelFixed <- function(model, data, iter, warmup, thin, init, seed, chain = 1,
                          stepsize = 1, adapt_delt = 0.8, max_depth = 10, refresh = 100){
  modelName <- basename(model)
  model <- file.path(model, modelName)
  print(paste0(model, " sample algorithm=fixed_param",
               " num_samples=1 num_warmup=0",
               " data file=", data,
               " random seed=", seed,
               " output file=", paste(model, chain, ".csv", sep = ""),
               " refresh=", refresh))
  
  system(paste0(model, " sample algorithm=fixed_param",
                " num_samples=1 num_warmup=0",
                " data file=", data,
                " init=", init,
                " random seed=", seed,
                " output file=", paste(model, chain, ".csv", sep = ""),
                " refresh=", refresh), invisible = FALSE)
  
}

################################################################################################
## Simulate ME-2 plasma concentrations and ANC values
## using mrgsolve.

nSub <- 15; # number of subjects
nIIV <- 7; # number of parameters with inter-individual variations

code <- '
$PARAM CL = 10, Q = 15, VC = 35, VP = 105, KA = 2.0, MTT = 125, 
       Circ0 = 5, alpha = 3E-4, gamma = 0.17, WT = 70

$SET delta=0.1 // simulation grid

$CMT GUT CENT PERI PROL TRANSIT1 TRANSIT2 TRANSIT3 CIRC

$MAIN
// Individual PK parameters
double CLi = exp(log(CL) + 0.75*log(WT/70) + ETA(1));
double Qi = exp(log(Q) + 0.75*log(WT/70) + ETA(2));
double VCi = exp(log(VC) + log(WT/70) + ETA(3));
double VPi = exp(log(VP) + log(WT/70) + ETA(4));

// Indivudal PD parameters
double MTTi = exp(log(MTT) + ETA(5));
double Circ0i = exp(log(Circ0) + ETA(6));
double alphai = exp(log(alpha) + ETA(7));

// Reparametrization
double k10 = CLi / VCi;
double k12 = Qi / VCi;
double k21 = Qi / VPi;
double ktr = 4/MTTi;

$ODE 
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (k10 + k12) * CENT + k21 * PERI;
dxdt_PERI = k12 * CENT - k21 * PERI;
dxdt_PROL = ktr * (PROL + Circ0i) * 
              ((1 - alphai * CENT/VCi) * pow(Circ0i/(CIRC + Circ0i),gamma) - 1);
dxdt_TRANSIT1 = ktr * (PROL - TRANSIT1);
dxdt_TRANSIT2 = ktr * (TRANSIT1 - TRANSIT2);
dxdt_TRANSIT3 = ktr * (TRANSIT2 - TRANSIT3);
dxdt_CIRC = ktr * (TRANSIT3 - CIRC);

$OMEGA name="IIV"
0.0025 0.0025 0.0025 0.0025 0.04 0.0256 0.016
// 0 0 0 0 0 0 0

$SIGMA 0.01 0.0025 // 0 0

$TABLE
double CP = CENT/VCi;
double DV1 = CENT/VCi * exp(EPS(1));
double DV2 = (CIRC + Circ0i) * exp(EPS(2));
double WEIGHT = WT;

$CAPTURE CP DV1 DV2 WEIGHT
'

mod <- mread("acum", tempdir(), code)
e1 <- expand.ev(amt = 80 * 1000, addl = 14, ii = 12, WT = rnorm(nSub, 70, 15))
out <- mod %>% data_set(e1) %>% carry.out(dose) %>% Req(CP,DV1,DV2) %>% mrgsim(end=500)
plot(out, DV1~time|factor(ID),scales="same")
plot(out, DV2~time|factor(ID),scales="same")

## Observation and dosing times
# doseTimes <- seq(0, 168, by = 12)
xpk <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
xpk <- c(xpk, xpk + 12, seq(24, 156, by = 12), c(xpk, 12, 18, 24) + 168)
xneut <- seq(0, 672, by = 24)
time <- sort(unique(c(xpk, xneut)))
# time <- sort(unique(c(xpk, xneut, doseTimes)))

## Assemble data set for Stan
xdata <- mod %>% data_set(e1) %>%
           carry.out(cmt, ii, addl, rate, amt, evid, ss) %>%
           mrgsim(Req = "DV1, DV2, WEIGHT", end = -1, add = time, rescort = 3) %>%
           as.data.frame

xdata <- xdata %>%
    mutate(DV1 = ifelse(time %in% xpk & time != 0 & evid == 0, DV1, NA),
           DV2 = ifelse(time %in% xneut & evid == 0, DV2, NA))

xdata$cmt[xdata$cmt == 0] <- 2 # switch from mrgsolve to NONMEM convention

nt <- nrow(xdata)

## Subject specific data
start <- (1:nt)[!duplicated(xdata$ID)]
end <- c(start[-1] -1, nt)
xsub <- subset(xdata, !duplicated(ID))
weight <- xsub$WEIGHT

## Check simulated data using plots
p1 <- ggplot(xdata %>% filter(!is.na(DV1)), aes(x = time, y = DV1))
p1 + geom_point() + geom_line() +
  labs(x = "time (h)", y = "ME-2 plasma concentration (ng/mL)") +
             theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                   legend.position = "none", strip.text = element_text(size = 8)) +
  facet_wrap(~ID)

p1 <- ggplot(xdata %>% filter(!is.na(DV2)), aes(x = time, y = DV2))
p1 + geom_point() + geom_line() +
    labs(x = "time (h)",
         y = "ANC") +
             theme(text = element_text(size = 12), axis.text = element_text(size = 12),
                   legend.position = "none", strip.text = element_text(size = 8)) +
  facet_wrap(~ID)

## Indices of records containing observed concentrations
iObsPK <- with(xdata, (1:nrow(xdata))[!is.na(DV1) & evid == 0])
nObsPK <- length(iObsPK)
## Indices of records containing observed neutrophil counts
iObsPD <- with(xdata, (1:nrow(xdata))[!is.na(DV2) & evid == 0])
nObsPD <- length(iObsPD)

## Creqte parameters for priors distribution
CLHatPrior = 10
QHatPrior = 15
V1HatPrior = 35
V2HatPrior = 105
kaHatPrior = 2
CLHatPriorCV = 0.10
QHatPriorCV = 0.18
V1HatPriorCV = 0.14
V2HatPriorCV = 0.17
kaHatPriorCV = 0.16

circ0HatPrior <- 5
circ0HatPriorCV <- 0.20
mttHatPrior <- 125
mttHatPriorCV <- 0.2
gammaPrior <- 0.17
gammaPriorCV <- 0.2
alphaHatPrior <- 2.0E-4
alphaHatPriorCV <- 0.2

## create data set
data <- with(xdata,
             list(
                 nt = nt,
                 nObsPK = nObsPK,
                 iObsPK = iObsPK,
                 nObsPD = nObsPD,
                 iObsPD = iObsPD,
                 amt = amt,
                 cmt = cmt,
                 evid = evid,
                 time = time,
                 ii = ii,
                 addl = addl,
                 ss = ss,
                 rate = rate,
                 cObs = DV1[iObsPK],
                 neutObs = DV2[iObsPD],
                 
                 nSubjects = nSub,
                 nIIV = nIIV,
                 start = start,
                 end = end,
                 weight = weight,
                 
                 CLHatPrior = CLHatPrior,
                 QHatPrior = QHatPrior,
                 V1HatPrior = V1HatPrior,
                 V2HatPrior = V2HatPrior,
                 kaHatPrior = kaHatPrior,
                 CLHatPriorCV = CLHatPriorCV,
                 QHatPriorCV = QHatPriorCV,
                 V1HatPriorCV = V1HatPriorCV,
                 V2HatPriorCV = V2HatPriorCV,
                 kaHatPriorCV = kaHatPriorCV,
                 circ0HatPrior = circ0HatPrior,
                 circ0HatPriorCV = circ0HatPriorCV,
                 mttHatPrior = mttHatPrior,
                 mttHatPriorCV = mttHatPriorCV,
                 gammaPrior = gammaPrior,
                 gammaPriorCV = gammaPriorCV,
                 alphaHatPrior = alphaHatPrior,
                 alphaHatPriorCV = alphaHatPriorCV
))

## create initial estimates
init <- function(){
    list(CLHat = exp(rnorm(1, log(CLHatPrior), CLHatPriorCV)),
         QHat = exp(rnorm(1, log(QHatPrior), QHatPriorCV)),
         V1Hat = exp(rnorm(1, log(V1HatPrior), V1HatPriorCV)),
         V2Hat = exp(rnorm(1, log(V2HatPrior), V2HatPriorCV)),
         kaHat = exp(rnorm(1, log(kaHatPrior), kaHatPriorCV)),
##         sigma = runif(1, 0.5, 2),
         sigma = 0.2,
         alphaHat = exp(rnorm(1, log(alphaHatPrior), alphaHatPriorCV)),
         mttHat = exp(rnorm(1, log(mttHatPrior), mttHatPriorCV)),
         circ0Hat = exp(rnorm(1, log(circ0HatPrior), circ0HatPriorCV)),
         gamma = exp(rnorm(1, log(gammaPrior), gammaPriorCV)),
##         sigmaNeut = runif(1, 0.5, 2))
         sigmaNeut = 0.2,
         L = diag(nIIV),
         omega = exp(rnorm(nIIV, log(0.05), 0.5)),
         # omega = c(0.05, 0.05, 0.05, 0.05, 0.2, 0.16, sqrt(0.016)),
         etaStd = matrix(rep(0, nIIV * nSub), nrow = nIIV))
}

## use true values for initial conditions
## Note: we cannot perform an exact deterministic test, because
## the inter-individual variations are random.
if (FALSE) {
 init <- function(){
   list(CLHat = 10,
        QHat = 15,
        V1Hat = 35,
        V2Hat = 105,
        kaHat = 2.0,
        # sigma = runif(1, 0.5, 2),
        # sigma = 2.5E-7,
        sigma = 0.1,
        alphaHat = 3E-4,
        mttHat = 125,
        circ0Hat = 5,
        gamma = 0.17,
        # sigmaNeut = runif(1, 0.5, 2))
        sigmaNeut = 0.05,
        # sigmaNeut = 2.5E-7,
        L = diag(nIIV),
        omega = c(0.05, 0.05, 0.05, 0.05, 0.2, 0.16, sqrt(0.016)),
        # omega = rep(1E-7, nIIV),
        etaStd = matrix(rep(0, nIIV * nSub), nrow = nIIV)) }
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "kaHat", "sigma",
                      "alphaHat", "mttHat", "circ0Hat", "gamma", "sigmaNeut",
                      "omega")

## Additional variables to monitor
otherRVs <- c("cObsCond", "neutObsCond") ##, "log_lik")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

################################################################################################
# run Stan
nChains <- 4
nPost <- 100 ## Number of post-burn-in samples per chain after thinning
nBurn <- 100 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

chains <- 1:nChains

if (FALSE) {
  ## Run model with fixed parameters 
  mclapply(chains,
           function(chain, model, data, iter, warmup, thin, init) {
             runModelFixed(model = model, data = data,
                           # iter = iter, warmup = warmup, 
                           thin = thin,
                           init = init, seed = sample(1:999999, 1),
                           chain = chain,
                           refresh = 1)
             },
           model = file.path(modelDir, modelName),
           data = file.path(tempDir, "data.R"),
           init = file.path(tempDir, "init.R"),
           iter = n.iter, warmup = n.burnin, thin = n.thin,
           mc.cores = min(nChains, detectCores()))
}

mclapply(chains,
         function(chain, model, data, iter, warmup, thin, init, tempDir = "temp"){
             tempDir <- file.path(tempDir, chain)
             dir.create(tempDir)
             with(data, stan_rdump(ls(data), file = file.path(tempDir, "data.R")))
             inits <- init()
             with(inits, stan_rdump(ls(inits), file = file.path(tempDir, "init.R")))
             runModel(model = model, data = file.path(tempDir, "data.R"),
                      iter = iter, warmup = warmup, thin = thin,
                      init = file.path(tempDir, "init.R"), seed = sample(1:999999, 1),
                      chain = chain, refresh = 1,
                      adapt_delta = 0.95, stepsize = 0.01)
            },
         model = file.path(modelDir, modelName),
         data = data,
         init = init,
         iter = nIter, warmup = nBurnin, thin = nThin,
         mc.cores = min(nChains, detectCores()))

fit <- read_stan_csv(file.path(modelDir, modelName, glue(modelName, chains, ".csv")))

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
##load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters
 
dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

mcmcHistory(fit, parametersToPlot)
mcmcDensity(fit, parametersToPlot, byChain = TRUE)
mcmcDensity(fit, parametersToPlot)
pairs(fit, pars = parametersToPlot)

ptable <- parameterTable(fit, parametersToPlot)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

if(any(grepl("log_lik", names(fit)))){
    logLik <- extract_log_lik(fit)
    waic1 <- waic(logLik)
    loo1 <- loo(logLik)
    waiclooTable <- data.frame(parameter = c("elpd_waic", "p_waic", "waic",
                                        "elpd_loo", "p_loo", "looic"),
                               estimate = c(with(waic1, c(elpd_waic, p_waic, waic)),
                                   with(loo1, c(elpd_loo, p_loo, looic))),
                               se = c(with(waic1, c(se_elpd_waic, se_p_waic, se_waic)),
                                   with(loo1, c(se_elpd_loo, se_p_loo, se_looic))))
    write.csv(waiclooTable, file.path(tabDir, paste(modelName, "WaicLoo.csv", sep = "")))
    print(waic1)
    print(loo1)
}

################################################################################################
## posterior predictive distributions

# prediction of future observations in the same studies, i.e., posterior predictions
# conditioned on observed data from the same study

## PK
pred <- as.data.frame(fit, pars = "cHatObsCond") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
                      median = quantile(value, probs = 0.5, na.rm = TRUE),
                      ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
                          bind_cols(xdata[!is.na(xdata$DV1),])

p1 <- ggplot(pred, aes(x = time, y = DV1))
    p1 <- p1 + geom_point() +
        labs(x = "time (h)",
             y = "ME-2 plasma concentration (ng/mL)") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12),
legend.position = "none", strip.text = element_text(size = 8)) + facet_wrap(~ID)

print(p1 + geom_line(aes(x = time, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))

## PD
pred <- as.data.frame(fit, pars = "neutHatObsCond") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
                      median = quantile(value, probs = 0.5, na.rm = TRUE),
                      ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
                          bind_cols(xdata[!is.na(xdata$DV2),])

p1 <- ggplot(pred, aes(x = time, y = DV2))
    p1 <- p1 + geom_point() +
        labs(x = "time (h)",
             y = "ANC") +
                 theme(text = element_text(size = 12), axis.text = element_text(size = 12),
legend.position = "none", strip.text = element_text(size = 8)) + facet_wrap(~ID)

print(p1 + geom_line(aes(x = time, y = median)) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))

## prediction of future observations in new subjects

## PK
pred <- as.data.frame(fit, pars = "cHatObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata[!is.na(xdata$DV1),])

p1 <- ggplot(pred, aes(x = time, y = DV1))
p1 <- p1 + geom_point() +
  labs(x = "time (h)",
       y = "ME-2 plasma concentration (ng/mL)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) + facet_wrap(~ID)

print(p1 + geom_line(aes(x = time, y = median)) +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))

## PD
pred <- as.data.frame(fit, pars = "neutHatObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata[!is.na(xdata$DV2),])

p1 <- ggplot(pred, aes(x = time, y = DV2))
p1 <- p1 + geom_point() +
  labs(x = "time (h)",
       y = "ANC") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) + facet_wrap(~ID)

print(p1 + geom_line(aes(x = time, y = median)) +
        geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))


dev.off()
