rm(list = ls())
gc()

modelName <- "ttpn2"
simModelName <- "ttpn2Sim"

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- scriptDir
figDir <- file.path(projectDir, "figure")
tabDir <- file.path(projectDir, "table")
dataDir <- file.path(projectDir, "derived")
modelDir <- file.path(projectDir)
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path(scriptDir, "../R/tools")

## .libPaths("lib")

library(rstan)
library(bayesplot)
## Go back to default ggplot2 theme that was overridden by bayesplot
theme_set(theme_gray())
library(tidyverse)
library(parallel)
library(sampling)
library(survival)
source(file.path(toolsDir, "stanTools.R"))
source(file.path(toolsDir, "functions.R"))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(11191951) ## not required but assures repeatable results

################################################################################################
### Simulate ME-2 plasma concentrations and ANC values

## Parameter values

CLHat = 0.0052
omegaCL = 0.25
VHat = 0.25
omegaV = 0.25
sigma = 0.15
omega = c(omegaCL, omegaV)
rho = diag(2)

ke0 = 3.60E-4
alpha = 2.26E-6
beta = 1.37

## Observation and dosing times
dt = 12
tObs <- seq(0, 6 * 21 * 24, by = dt)

## Assemble data set for simulation using Stan
obsData <- data.frame(time = tObs) %>%
    mutate(cmt = 1,
           evid = 0,
           addl = 0,
           ii = 0)

doseData <- data.frame(time = 0) %>%
    mutate(cmt = 1,
           evid = 1,
           addl = 5,
           ii = 21 * 24)

nId = 20 * 3
allData <- doseData %>%
    bind_rows(obsData) %>%
    merge(data.frame(id = 1:nId, amt = rep(c(1.2, 1.8, 2.4), ea = nId / 3))) %>%
    arrange(id, time, desc(evid)) %>%
    mutate(amt = ifelse(evid == 0, 0, amt),
           rate = amt)

nt <- nrow(allData)
start <- (1:nt)[!duplicated(allData$id)]
end <- c(start[-1] - 1, nt)

## Generate sparse PK sampling times by randomly selecting 6 sample times per patient
iPKObs = sampling:::strata(allData, stratanames = "id", method = "srswor",
                size = rep(6, length(unique(allData$id))))$ID_unit
nPKObs = length(iPKObs)

dataSim <- with(allData,
                list(nId = nId,
                     nt = nt,
                     nPKObs = nPKObs,
                     iPKObs = iPKObs,
                     amt = 1000 * amt,
                     rate = 0 * rate,
                     cmt = cmt,
                     evid = evid,
                     ii = ii,
                     addl = addl,
                     time = time,
                     start = start,
                     end = end,
                     CLHat = CLHat,
                     VHat = VHat,
                     ke0 = ke0,
                     alpha = alpha,
                     beta = beta,
                     nRandom = length(omega),
                     omega = omega,
                     rho = rho,
                     sigma = sigma))

parameters = c("cObs", "cHat", "ce", "cdf", "CL", "V", "x")

### Simulate using Stan

sim <- stan(file = file.path(modelDir, paste(simModelName, ".stan", sep = "")),
            data = dataSim,
            pars = parameters,
            algorithm = "Fixed_param",
            iter = 1,
            chains = 1)

################################################################################################
### Assemble data set for fitting via Stan

dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

pkData <- allData %>%
  bind_cols(as.data.frame(sim, pars = "cHat") %>%
            gather(factor_key = TRUE) %>%
            select(cObs = value)) %>%
  bind_cols(as.data.frame(sim, pars = "ce") %>%
            gather(factor_key = TRUE) %>%
            select(ce = value)) %>%
  mutate(type = "pk",
         censor = FALSE)

pkPar <- as.data.frame(sim, pars = "CL")  %>%
    gather(factor_key = TRUE) %>%
    select(CL = value) %>%
    bind_cols(as.data.frame(sim, pars = "V")  %>%
    gather(factor_key = TRUE) %>%
    select(V = value)) %>%
    mutate(id = 1:nId)

cdf <- allData[c("id", "time")] %>%
  bind_cols(as.data.frame(sim, pars = "cdf") %>%
              gather(factor_key = TRUE) %>%
              select(cdf = value))

## Simulate time to PN using simulated cdf.
tpn <- cdf %>%
  group_by(id) %>%
  summarize(tpn = approx(cdf, time, xout = runif(1))$y,
            censor = is.na(tpn),
            tpn = ifelse(censor, max(time), tpn))

tpn = pkData %>%
    filter(!duplicated(id)) %>%
    select(-time, -cObs, -censor) %>%
    left_join(tpn) %>%
    rename(time = tpn) %>%
    mutate(cObs = NA,
           type = "tpn",
           evid = 0,
           addl = 0,
           ii = 0,
           amt = 0,
           rate = 0)

doseData <- allData %>%
    filter(evid == 1) %>%
    mutate(cObs = NA,
           censor = FALSE,
           type = "dose")

xdata <- pkData %>%
    bind_rows(tpn) %>%
    bind_rows(doseData) %>%
    left_join(pkPar) %>%
    arrange(id, time, desc(evid))

doseData <- doseData %>%
    select(id, dose = amt)

xdata <- xdata %>%
    left_join(doseData)

doses <- sort(unique(xdata$dose))
for(thisDose in doses){
  xplot <- xdata %>% filter(dose == thisDose & type == "pk")
  p1 <- ggplot(xplot %>% filter(!is.na(cObs)), aes(x = time, y = cObs))
  p1 <- p1 + geom_line() +
    labs(title = paste("dose =", thisDose, "mg"),
         x = "time (h)",
         y = "plasma drug concentration (ng/mL)") +
    theme(text = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1)

  p1 <- ggplot(xplot %>% filter(!is.na(ce)), aes(x = time, y = ce))
  p1 <- p1 + geom_line() +
    labs(title = paste("dose =", thisDose, "mg"),
         x = "time (h)",
         y = "effect site concentration (ng/mL)") +
    theme(text = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1)
}

tpn <- tpn %>% left_join(doseData)
pnSurv <- with(tpn, Surv(time = time, event = !censor))
fitTpn <- survfit(pnSurv ~ -1 + dose, tpn, conf.int = 0.9)
plot(fitTpn, fun = "event", ylim = c(0, 1), col = 1:3)
legend(100, .9, c("1.2 mg/kg", "1.8 mg/kg", "2.4 mg/kg"), lty = 1, col = 1:3)

xdata <- xdata %>%
    filter(type %in% c("dose", "tpn"))

dir.create(dataDir)
write.csv(xdata, file = file.path(dataDir, paste(modelName, "Data.csv", sep = "")))

## Indices of records containing observed time to PN
iPNObs <- with(xdata, (1:nrow(xdata))[type == "tpn" & !censor])
nPNObs <- length(iPNObs)
## Indices of records containing censored time to PN
iPNCens <- with(xdata, (1:nrow(xdata))[type == "tpn" & censor])
nPNCens <- length(iPNCens)

nt <- nrow(xdata)
start <- (1:nt)[!duplicated(xdata$id)]
end <- c(start[-1] - 1, nt)

data <- with(xdata,
                list(nId = nId,
                     nt = nt,
                     nPNObs = nPNObs,
                     iPNObs = iPNObs,
                     nPNCens = nPNCens,
                     iPNCens = iPNCens,
                     amt = 1000 * amt,
                     rate = 0 * rate,
                     cmt = cmt,
                     evid = evid,
                     ii = ii,
                     addl = addl,
                     time = time,
                     start = start,
                     end = end,
                     CL = CL[!duplicated(id)],
                     V = V[!duplicated(id)]
                     ))

### create initial estimates
init <- function(){
  list(alpha = exp(rnorm(1, log(2.0E-6), 0.5)),
       beta = 1 + exp(rnorm(1, log(0.5), 0.5)),
       ke0 = exp(rnorm(1, log(4.0E-4), 0.5))
  )
}

### Specify the variables for which you want history and density plots

parametersToPlot <- c("ke0", "alpha", "beta")

## Additional variables to monitor
otherRVs <- c("cdfPred")

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
            cores = nChains,
            refresh = 10,
            control = list(adapt_delta = 0.95, stepsize = 0.01))

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
##load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters

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

ntPred <- 253
dt <- 12
tPred = dt * 0:(ntPred - 1)
nSamp <- nChains * nPost

cdfPred <- as.matrix(fit, pars = "cdfPred")
dim(cdfPred) <- c(nSamp, nId, ntPred)

dose1 <- cdfPred[,1:(nId / 3),]
dim(dose1) <- c(nSamp * (nId / 3), ntPred)
dose2 <- cdfPred[,nId / 3 + 1:(nId / 3),]
dim(dose2) <- c(nSamp * (nId / 3), ntPred)
dose3 <- cdfPred[,2 * nId / 3 + 1:(nId / 3),]
dim(dose3) <- c(nSamp * (nId / 3), ntPred)

pred <- as.data.frame(dose1) %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  mutate(dose = doses[1],
         time = tPred) %>%
  bind_rows(
    as.data.frame(dose2) %>%
      gather(factor_key = TRUE) %>%
      group_by(key) %>%
      summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
                median = quantile(value, probs = 0.5, na.rm = TRUE),
                ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
      mutate(dose = doses[2],
             time = tPred)
  ) %>%
  bind_rows(
    as.data.frame(dose3) %>%
      gather(factor_key = TRUE) %>%
      group_by(key) %>%
      summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
                median = quantile(value, probs = 0.5, na.rm = TRUE),
                ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
      mutate(dose = doses[3],
             time = tPred)
  )

tpnObs <- data.frame(time = fitTpn$time,
                     ecdf = 1 - fitTpn$surv,
                     dose = rep(doses, fitTpn$strata))


p1 <- ggplot(pred, aes(x = time, y = median))
p1 <- p1 + geom_line() +
  labs(x = "time (h)",
       y = "Pr(PN grade 2+)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) +
  facet_wrap(~ dose)
p1 <- p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)
p1 + geom_step(data = tpnObs, aes(x = time, y = ecdf))

dev.off()

