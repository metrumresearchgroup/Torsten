rm(list = ls())
gc()

setwd("~/Desktop/Code/torsten/example-models/PKPD/torsten/R")
.libPaths("~/svn-StanPmetrics/script/lib")
modelName <- "effCpt"

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

# source(file.path(scriptDir, "pkgSetup.R"))
.libPaths("~/svn-StanPmetrics/script/lib")

# library(rstan)
# library(metrumrg)
# library(ggplot2)
# library(plyr)
# library(dplyr)
# library(tidyr)
# library(loo)
# library(parallel)
# source(file.path(toolsDir, "stanTools.R"))
# source(file.path(toolsDir, "functions.R"))

rstan_options(auto_write = TRUE)

set.seed(10271998) ## not required but assures repeatable results

################################################################################################

## get data file
library(dplyr)
xdata1 <- read.csv(file.path(dataDir, "phase1effcpt.csv"), as.is = TRUE)
xdata2 <- read.csv(file.path(dataDir, "phase2effcpt.csv"), as.is = TRUE)

xdata1 <- xdata1 %>%
  select(id = subject, time, weight, dose, cobs, resp = fxa.inh.obs) %>%
  mutate(cobs = as.numeric(cobs),
         study = 1,
         evid = 0) %>%
  filter(dose > 0, time > 0)

dose1 <- xdata1 %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(time = 0,
         evid = 1,
         cobs = NA,
         resp = NA)

xdata1 <- xdata1 %>%
  bind_rows(dose1) %>%
  arrange(id, time, desc(evid))

xdata2 <- xdata2 %>%
  filter(drug == "ME-2") %>%
  select(id = patient, time, weight, dose, cobs, resp = fxa.inh.obs) %>%
  mutate(cobs = as.numeric(cobs),
         study = 2,
         evid = 0,
         id = id + max(xdata1$id)) %>%
  filter(time > 0)

dose2 <- xdata2 %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(evid = 1,
         cobs = NA,
         resp = NA) %>%
  select(-time) %>%
  merge(data.frame(time = seq(0, 13 * 12, by = 12)))

xdata2 <- xdata2 %>%
  bind_rows(dose2) %>%
  arrange(id, time, desc(evid))

xdata <- xdata1 %>%
  bind_rows(xdata2)

nt <- nrow(xdata)
start <- (1:nt)[!duplicated(xdata$id)]
end <- c(start[-1] - 1, nt)

## Indices of records containing observed concentrations
iObs <- with(xdata, (1:nrow(xdata))[!is.na(cobs) & evid == 0])
nObs <- length(iObs)

## create data set
data <- with(xdata,
             list(
               nSubjects = length(unique(id)),
               nt = nt,
               nObs = nObs,
               iObs = iObs,
               amt = dose,
               cmt = rep(1, nt),
               evid = evid,
               start = start,
               end = end,
               time = time,
               cObs = cobs[iObs],
               respObs = resp[iObs],
               weight = weight[!duplicated(id)],
               rate = rep(0, nt),
               ii = rep(0, nt),
               addl = rep(0, nt),
               ss = rep(0, nt)
             ))

## create initial estimates
init <- function(){
  list(CLHat = exp(rnorm(1, log(10), 0.2)),
       QHat = exp(rnorm(1, log(20), 0.2)),
       V1Hat = exp(rnorm(1, log(70), 0.2)),
       V2Hat = exp(rnorm(1, log(70), 0.2)),
       kaHat = exp(rnorm(1, log(2), 0.2)),
       ke0Hat = exp(rnorm(1,log(1),0.2)),
       EC50Hat = exp(rnorm(1,log(100),0.2)),
       omega = exp(rnorm(5, log(0.25), 0.5)),
       rho = diag(5),
       omegaKe0 = exp(rnorm(1, log(0.25), 0.5)),
       omegaEC50 = exp(rnorm(1, log(0.25), 0.5)),
       sigma = 0.5,
       sigmaResp = 20,
       logtheta = matrix(rep(log(c(10, 20, 70, 70, 2)), ea = 200), nrow = 200),
       logKe0 = rep(log(1), 200),
       logEC50 = rep(log(100), 200))
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "kaHat",
                      "ke0Hat", "EC50Hat", 
                      "sigma", "sigmaResp", "omega", "rho",
                      "omegaKe0", "omegaEC50")

## Additional variables to monitor
otherRVs <- c("cObsCond", "cObsPred", "respObsCond", "respObsPred",
              "CL", "Q", "V1", "V2", "ka", "ke0", "EC50") ##, "log_lik")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

################################################################################################
# run Stan

nChains <- 4
nPost <- 1000 ## Number of post-burn-in samples per chain after thinning
nBurn <- 1000 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- nPost * nThin
nBurnin <- nBurn * nThin

RNGkind("L'Ecuyer-CMRG")
library(parallel)
mc.reset.stream()

source(file.path(toolsDir, "cmdStanTools.R"))
source(file.path(toolsDir, "stanTools.R"))
compileModel(model = file.path(modelDir, modelName), stanDir = stanDir)

chains <- 1:nChains
library(rstan)
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


library(metrumrg)
fit <- read_stan_csv(file.path(modelDir, modelName, glue(modelName, chains, ".csv")))

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))
load(file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

################################################################################################
## posterior distributions of parameters

dir.create(figDir)
dir.create(tabDir)

## open graphics device
pdf(file = file.path(figDir, paste(modelName,"Plots%03d.pdf", sep = "")),
    width = 6, height = 6, onefile = F)

## Remove diagonal & redundant elements of rho
dimRho <- nrow(init()$rho)
parametersToPlot <- c(parametersToPlot,
                      paste("rho[", matrix(apply(expand.grid(1:dimRho, 1:dimRho), 1, paste, collapse = ","),
                                           ncol = dimRho)[upper.tri(diag(dimRho), diag = FALSE)], "]", sep = ""))
parametersToPlot <- setdiff(parametersToPlot, "rho")

library(tidyr)
mcmcHistory(fit, parametersToPlot)
mcmcDensity(fit, parametersToPlot, byChain = TRUE)
mcmcDensity(fit, parametersToPlot)
pairs(fit, pars = parametersToPlot[!grepl("rho", parametersToPlot)])

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

pred <- as.data.frame(fit, pars = "cObsCond") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

xplot <- subset(pred, study == 1)
doses <- sort(unique(xplot$dose))
for(thisDose in doses){
  p1 <- ggplot(subset(xplot, dose == thisDose), aes(x = time, y = cobs))
  p1 <- p1 + geom_point() +
    labs(title = paste("study 1", thisDose, "mg\n individual predictions"),
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

xplot <- subset(pred, study == 2)

nPerPage = 25
subjects <- sort(unique(xplot$id))
nSubjects <- length(subjects)
nPages <- ceiling(nSubjects / nPerPage)
subjects <- data.frame(id = subjects,
                       page = sort(rep(1:nPages, length = nSubjects)),
                       stringsAsFactors = FALSE)
xplot <- stableMerge(xplot, subjects)

for(thisPage in 1:nPages){
  p1 <- ggplot(subset(xplot, page == thisPage), aes(x = time, y = cobs))
  p1 <- p1 + geom_point() +
    labs(title = "study 2 20 mg\n individual predictions",
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
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

xplot <- subset(pred, study == 1)
doses <- sort(unique(xplot$dose))
for(thisDose in doses){
  p1 <- ggplot(subset(xplot, dose == thisDose), aes(x = time, y = resp))
  p1 <- p1 + geom_point() +
    labs(title = paste("study 1", thisDose, "mg\n individual predictions"),
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

xplot <- subset(pred, study == 2)

nPerPage = 25
subjects <- sort(unique(xplot$id))
nSubjects <- length(subjects)
nPages <- ceiling(nSubjects / nPerPage)
subjects <- data.frame(id = subjects,
                       page = sort(rep(1:nPages, length = nSubjects)),
                       stringsAsFactors = FALSE)
xplot <- stableMerge(xplot, subjects)

for(thisPage in 1:nPages){
  p1 <- ggplot(subset(xplot, page == thisPage), aes(x = time, y = resp))
  p1 <- p1 + geom_point() +
    labs(title = "study 2 20 mg\n individual predictions",
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

# prediction of future observations in a new study, i.e., posterior predictive distributions

pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05, na.rm = TRUE),
            median = quantile(value, probs = 0.5, na.rm = TRUE),
            ub = quantile(value, probs = 0.95, na.rm = TRUE)) %>%
  bind_cols(xdata)

xplot <- subset(pred, study == 1)
doses <- sort(unique(xplot$dose))
for(thisDose in doses){
  p1 <- ggplot(subset(xplot, dose == thisDose), aes(x = time, y = cobs))
  p1 <- p1 + geom_point() +
    labs(title = paste("study 1", thisDose, "mg\n population predictions"),
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

xplot <- subset(pred, study == 2)

nPerPage = 25
subjects <- sort(unique(xplot$id))
nSubjects <- length(subjects)
nPages <- ceiling(nSubjects / nPerPage)
subjects <- data.frame(id = subjects,
                       page = sort(rep(1:nPages, length = nSubjects)),
                       stringsAsFactors = FALSE)
xplot <- stableMerge(xplot, subjects)

for(thisPage in 1:nPages){
  p1 <- ggplot(subset(xplot, page == thisPage), aes(x = time, y = cobs))
  p1 <- p1 + geom_point() +
    labs(title = "study 2 20 mg\n population predictions",
         x = "time (h)",
         y = "ME-2 plasma concentration (ng/mL)") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
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

xplot <- subset(pred, study == 1)
doses <- sort(unique(xplot$dose))
for(thisDose in doses){
  p1 <- ggplot(subset(xplot, dose == thisDose), aes(x = time, y = resp))
  p1 <- p1 + geom_point() +
    labs(title = paste("study 1", thisDose, "mg\n population predictions"),
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

xplot <- subset(pred, study == 2)

nPerPage = 25
subjects <- sort(unique(xplot$id))
nSubjects <- length(subjects)
nPages <- ceiling(nSubjects / nPerPage)
subjects <- data.frame(id = subjects,
                       page = sort(rep(1:nPages, length = nSubjects)),
                       stringsAsFactors = FALSE)
xplot <- stableMerge(xplot, subjects)

for(thisPage in 1:nPages){
  p1 <- ggplot(subset(xplot, page == thisPage), aes(x = time, y = resp))
  p1 <- p1 + geom_point() +
    labs(title = "study 2 20 mg\n population predictions",
         x = "time (h)",
         y = "response") +
    theme(text = element_text(size = 12), axis.text = element_text(size = 12),
          legend.position = "none", strip.text = element_text(size = 8)) +
    facet_wrap(~ id)
  print(p1 + geom_line(aes(x = time, y = median)) +
          geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25))
}

dev.off()
