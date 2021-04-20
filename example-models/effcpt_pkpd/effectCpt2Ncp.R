rm(list = ls())
gc()

modelName <- "effectCpt2Ncp"
simModelName <- "effectCpt2NcpSim"

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- scriptDir
figDir <- file.path(projectDir, "deliv", "figure", modelName)
tabDir <- file.path(projectDir, "deliv", "table", modelName)
dataDir <- file.path(projectDir, "data", "derived")
modelDir <- projectDir
outDir <- file.path(modelDir, modelName)
## toolsDir <- file.path(scriptDir, "tools")

.libPaths("~/.R/R_libs")

## tools
qnorm.trunc = function(p,mean=0,sd=1,lower=-Inf,upper=Inf)
	qnorm(p*pnorm(upper,mean,sd)+(1-p)*pnorm(lower,mean,sd),mean,sd)

rnorm.trunc = function(n,mean=0,sd=1,lower=-Inf,upper=Inf)
	qnorm.trunc(runif(n),mean,sd,lower,upper)

library(cmdstanr)
set_cmdstan_path("../../cmdstan")
library(bayesplot)
## Go back to default ggplot2 theme that was overridden by bayesplot
## theme_set(theme_gray())
library(tidyverse)
## library(parallel)
## source(file.path(toolsDir, "stanTools.R"))
## source(file.path(toolsDir, "functions.R"))

set.seed(11191951) ## not required but assures repeatable results

################################################################################################
### Simulate ME-2 plasma concentrations and ANC values

## Parameter values
ka <- 1.5
CL <- 10 # L/h
V <- 35 # L
sigma <- 0.15

ke0 <- 0.5
E0 <- 80
Emax <- 40
EC50 <- 250
gamma <- 2

sigmaPD <- 0.05

omega <- c(0.25, 0.4, 0.25, 0.25, 0.25, 0.25)
rho <- diag(6)

## Observation and dosing times
doses = c(10, 20, 40, 80)
xpk <- c(0, 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2,3,4,6,8)
## xpk <- c(xpk, xpk + 12, seq(24, 156, by = 12), c(xpk, 12, 18, 24) + 168)
xpk <- c(xpk, xpk + 12, seq(24, 120, by = 12))
xpd <- xpk[xpk %% 24 == 0]
time <- sort(unique(c(xpk, xpd)))

nPerDose <- 3
nId <- nPerDose * length(doses) ## Number of individuals
weight = rnorm.trunc(nId, 70, 15, 50, 100)

## Assemble data set for simulation using Stan
obsData <- data.frame(id = 1:nId,
                      amt = rep(doses * 1000, nPerDose)) %>%
  merge(data.frame(time = time)) %>%
  mutate(rate = 0,
         cmt = 2,
         evid = 0,
         ii = 0,
         addl = 0,
         ss = 0)

doseData <- data.frame(id = 1:nId,
                       amt = rep(doses * 1000, nPerDose)) %>%
    mutate(time = 0,
           rate = 0,
           cmt = 1,
           evid = 1,
           ii = 12,
           addl = 14,
           ss = 0)

allData <- doseData %>%
  bind_rows(obsData) %>%
    arrange(id, time, desc(evid))

nt <- nrow(allData)
start <- (1:nt)[!duplicated(allData$id)]
end <- c(start[-1] - 1, nt)

dataSim <- with(allData,
                list(nId = nId,
                     nt = nt,
                     amt = amt,
                     rate = rate,
                     cmt = cmt,
                     evid = evid,
                     ii = ii,
                     addl = addl,
                     rate = rate,
                     ss = ss,
                     time = time,
                     start = start,
                     end = end,
                     weight = weight,
                     CLHat = CL,
                     VHat = V,
                     kaHat = ka,
                     ke0Hat = ke0,
                     EmaxHat = Emax,
                     EC50 = EC50,
                     E0Hat = E0,
                     gamma = gamma,
                     nRandom = length(omega),
                     omega = omega,
                     rho = rho,
                     sigma = sigma,
                     sigmaPD = sigmaPD))

### Simulate using Stan
mod.sim <- cmdstan_model(file.path(modelDir, paste(simModelName, ".stan", sep = "")),force_recompile=TRUE,quiet=FALSE)
sim <- mod.sim$sample(data=dataSim, fixed_param=TRUE, seed=3829, iter_sampling = 1)

################################################################################################
### Assemble data set for fitting via Stan

xdata <- allData %>%
    bind_cols(as.data.frame(sim$draws(variables="cObs")) %>%
              gather(factor_key = TRUE) %>%
              select(cObs = value)) %>%
    bind_cols(as.data.frame(sim$draws(variables = "respObs")) %>%
              gather(factor_key = TRUE) %>%
              select(respObs = value))

xdata <- xdata %>%
    mutate(cObs = ifelse(time %in% xpk & time != 0 & evid == 0, cObs, NA),
           respObs = ifelse(time %in% xpd & evid == 0, respObs, NA))

head(xdata)

dir.create(figDir, recursive=TRUE)
dir.create(tabDir, recursive=TRUE)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
	width = 6, height = 6, onefile = F)

for(thisDose in doses){
  xplot <- xdata %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot %>% filter(!is.na(cObs)), aes(x = time, y = cObs))
  p1 <- p1 + geom_point() + geom_line() +
    labs(title = paste("dose =", thisDose, "mg"),
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    expand_limits(y = 0) +
    facet_wrap(~ id)
  print(p1)

  p1 <- ggplot(xplot %>% filter(!is.na(respObs)),
               aes(x = time, y = respObs))
  p1 <- p1 + geom_point() + geom_line() +
    labs(title = paste("dose =", thisDose, "mg"),
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    expand_limits(y = 0) +
    facet_wrap(~ id)
  print(p1)
}

## Indices of records containing observed concentrations
iObsPK <- with(xdata, (1:nrow(xdata))[!is.na(cObs) & evid == 0])
nObsPK <- length(iObsPK)
## Indices of records containing observed response
iObsPD <- with(xdata, (1:nrow(xdata))[!is.na(respObs) & evid == 0])
nObsPD <- length(iObsPD)

## Parameters for informative priors

CLPrior <- 10
VPrior <- 35
kaPrior <- 2
CLPriorCV <- 1 ## 0.10
VPriorCV <- 1 ## 0.14
kaPriorCV <- 1 ## 0.16

ke0Prior <- 0.5
ke0PriorCV <- 1
E0Prior <- 80
E0PriorCV <- 1
EmaxPrior <- 40
EmaxPriorCV <- 1
EC50Prior <- 250
EC50PriorCV <- 1
gammaPrior <- 2
gammaPriorCV <- 1

## create data set
data <- with(xdata,
             list(nId = nId,
                  nt = nt,
                  nObsPK = nObsPK,
                  iObsPK = iObsPK,
                  nObsPD = nObsPD,
                  iObsPD = iObsPD,
                  amt = amt,
                  cmt = cmt,
                  evid = evid,
                  ii = ii,
                  addl = addl,
                  rate = rate,
                  ss = ss,
                  time = time,
                  start = start,
                  end = end,
                  weight = weight,
                  cObs = cObs[iObsPK],
                  respObs = respObs[iObsPD],
                  CLPrior = CLPrior,
                  VPrior = VPrior,
                  kaPrior = kaPrior,
                  CLPriorCV = CLPriorCV,
                  VPriorCV = VPriorCV,
                  kaPriorCV = kaPriorCV,
                  ke0Prior = ke0Prior,
                  ke0PriorCV = ke0PriorCV,
                  E0Prior = E0Prior,
                  E0PriorCV = E0PriorCV,
                  EmaxPrior = EmaxPrior,
                  EmaxPriorCV = EmaxPriorCV,
                  EC50Prior = EC50Prior,
                  EC50PriorCV = EC50PriorCV,
                  gammaPrior = gammaPrior,
                  gammaPriorCV = gammaPriorCV
             ))
with(data, rstan::stan_rdump(ls(data), file = paste0(modelName,".data.R")))

### create initial estimates
init <- function(){
    CLHat = exp(rnorm(1, log(CLPrior), CLPriorCV))
    VHat = exp(rnorm(1, log(VPrior), VPriorCV))
    list(CLHat = CLHat,
         VHat = VHat,
         kaHat = CLHat / VHat + exp(rnorm(1, log(kaPrior), kaPriorCV)),
         sigma = 0.2,
         ke0Hat = exp(rnorm(1, log(ke0Prior), ke0PriorCV)),
         E0Hat = exp(rnorm(1, log(E0Prior), E0PriorCV)),
         EmaxHat = exp(rnorm(1, log(EmaxPrior), EmaxPriorCV)),
         EC50 = exp(rnorm(1, log(EC50Prior), EC50PriorCV)),
         gamma = exp(rnorm(1, log(gammaPrior), gammaPriorCV)),
         sigmaPD = 0.2,
         omega = exp(rnorm(6, rep(log(0.25), 6), 0.5)),
         L = diag(6),
         etaStd = matrix(rep(0, 6 * nId), nrow = 6)
    )
}

for (i in 1:4) {
    inits <- init()
    with(inits, rstan::stan_rdump(ls(inits), file = paste0("init.",i, ".R")))
}


### Specify the variables for which you want history and density plots

parametersToPlot <- c("CLHat", "VHat", "kaHat",
                      "sigma", "ke0Hat", "E0Hat", "EmaxHat", "EC50",
                      "gamma", "sigmaPD", "omega", "rho")

## Additional variables to monitor
otherRVs <- c("cObsCond", "respObsCond", "cObsPred", "respObsPred",
              "theta")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

nChains <- 4
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste(modelName, ".stan", sep = "")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin,
            init = init,
            chains = nChains,
            refresh = 10,
            control = list(adapt_delta = 0.9,
                           max_treedepth = 20))

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
##load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters

## Remove diagonal & redundant elements of rho
dimRho <- nrow(init()$L)
parametersToPlot <- c(parametersToPlot,
                     paste("rho[", matrix(apply(expand.grid(1:dimRho, 1:dimRho), 1, paste, collapse = ","),
                                          ncol = dimRho)[upper.tri(diag(dimRho), diag = FALSE)], "]", sep = ""))
parametersToPlot <- setdiff(parametersToPlot, "rho")

options(bayesplot.base_size = 12,
        bayesplot.base_family = "sans")
color_scheme_set(scheme = "brightblue")
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12))

rhats <- rhat(fit, pars = parametersToPlot)
mcmc_rhat(rhats) + yaxis_text() + myTheme

ratios1 <- neff_ratio(fit, pars = parametersToPlot)
mcmc_neff(ratios1) + yaxis_text() + myTheme

mcmcHistory(fit, pars = parametersToPlot, nParPerPage = 5, myTheme = myTheme)
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16, byChain = TRUE,
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))
mcmcDensity(fit, pars = parametersToPlot, nParPerPage = 16,
            myTheme = theme(text = element_text(size = 12), axis.text = element_text(size = 10)))

pairs(fit, pars = parametersToPlot[!grepl("rho", parametersToPlot)])

ptable <- monitor(as.array(fit, pars = parametersToPlot), warmup = 0, print = FALSE)
write.csv(ptable, file = file.path(tabDir, paste(modelName, "ParameterTable.csv", sep = "")))

################################################################################################
### posterior predictive distributions

#### prediction of future observations in the same studies, i.e., posterior predictions conditioned on observed data from the same study

pred <- as.data.frame(fit, pars = "cObsCond") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
                      median = quantile(value, probs = 0.5, na.rm = TRUE),
                      ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
                          bind_cols(xdata)

for(thisDose in doses){
  xplot <- pred %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot, aes(x = time, y = cObs))
  p1 <- p1 + geom_point() +
    labs(title = paste("individual predictions\n", "dose =", thisDose, "mg"),
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    expand_limits(y = 0) +
    facet_wrap(~ id)

  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

pred <- as.data.frame(fit, pars = "respObsCond") %>%
    gather(factor_key = TRUE) %>%
        group_by(key) %>%
            summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
                      median = quantile(value, probs = 0.5, na.rm = TRUE),
                      ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
                          bind_cols(xdata)

for(thisDose in doses){
  xplot <- pred %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot, aes(x = time, y = respObs))
  p1 <- p1 + geom_point() +
    labs(title = paste("individual predictions\n", "dose =", thisDose, "mg"),
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    expand_limits(y = 0) +
    facet_wrap(~ id)

  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

################################################################################################
### posterior predictive distributions

#### prediction of future observations in a new study, i.e., posterior predictive distributions

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

for(thisDose in doses){
  xplot <- pred %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot, aes(x = time, y = cObs))
  p1 <- p1 + geom_point() +
    labs(title = paste("population predictions\n", "dose =", thisDose, "mg"),
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    expand_limits(y = 0) +
    facet_wrap(~ id)

  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

pred <- as.data.frame(fit, pars = "respObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

for(thisDose in doses){
  xplot <- pred %>% filter(amt / 1000 == thisDose)
  p1 <- ggplot(xplot, aes(x = time, y = respObs))
  p1 <- p1 + geom_point() +
    labs(title = paste("population predictions\n", "dose =", thisDose, "mg"),
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    expand_limits(y = 0) +
    facet_wrap(~ id)

  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

dev.off()
