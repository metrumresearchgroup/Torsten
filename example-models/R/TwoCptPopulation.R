rm(list = ls())
gc()

modelName <- "TwoCptModelPopulation"

## Adjust directories to your settings.
scriptDir <- getwd()
projectDir <- dirname(scriptDir)
figDir <- file.path("deliv", "figure", modelName)
tabDir <- file.path("deliv", "table", modelName)
modelDir <- file.path(projectDir, modelName)
outDir <- file.path(modelDir, modelName)
toolsDir <- file.path("tools")

library(rstan)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)

source(file.path(toolsDir, "stanTools.R"))

rstan_options(auto_write = TRUE)
set.seed(11191951) ## not required but assures repeatable results

## read data
data <- read_rdump(file.path(modelDir, paste0(modelName,".data.R")))

## create initial estimates
nIIV <- 5
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

## Specify the variables for which you want history and density plots
parametersToPlot <- c("CLHat", "QHat", "V1Hat", "V2Hat", "kaHat", "sigma", "omega")

## Additional variables to monitor
otherRVs <- c("cObsCond", "cObsPred", "thetaM")

parameters <- c(parametersToPlot, otherRVs)
parametersToPlot <- c("lp__", parametersToPlot)

################################################################################################
## run Stan
nChains <- 1
nPost <- 1 ## Number of post-burn-in samples per chain after thinning
nBurn <- 0 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- nPost* nThin
nBurnin <- nBurn * nThin

fit <- stan(file = file.path(modelDir, paste0(modelName, ".stan")),
            data = data,
            pars = parameters,
            iter = nIter,
            warmup = nBurnin,
            thin = nThin,
            init = init,
            chains = nChains,
            cores = min(nChains, parallel::detectCores()))

dir.create(outDir)
save(fit, file = file.path(outDir, paste(modelName, "Fit.Rsave", sep = "")))

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

################################################################################################
## posterior predictive plots
data <- read_rdump(file.path(modelDir, paste0(modelName,".data.R")))
ID <- c()
for(i in 1:data$nSubjects) {
  if(i == data$nSubjects) data$start[i+1] <- length(data$time) + 1
  newID <- rep(i, data$start[i+1] - data$start[i] - 1)
  ID <- c(ID, newID)
}
data <- data.frame(data$cObs, data$time[(data$evid == 0)], ID)
data <- plyr::rename(data, c("data.cObs" = "cObs", "data.time..data.evid....0.." = "time"))

## predictions for furture observations in the same patient
pred <- as.data.frame(fit, pars = "cObsCond") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(title = "Individual Predictions", x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) +
  facet_wrap(~ ID)


## predictions for future observations in new subjects
pred <- as.data.frame(fit, pars = "cObsPred") %>%
  gather(factor_key = TRUE) %>%
  group_by(key) %>%
  summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data)

p1 <- ggplot(pred, aes(x = time, y = cObs))
p1 <- p1 + geom_point() +
  labs(title = "Population Predictions", x = "time (h)", y = "plasma concentration (mg/L)") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 8)) 
p1 + geom_line(aes(x = time, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25) +
  facet_wrap(~ ID)

dev.off()
