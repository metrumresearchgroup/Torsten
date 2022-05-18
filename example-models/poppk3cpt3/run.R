modelName <- "pk3cpt3"
scriptName <- paste(modelName, "Rmd", sep = ".")
fitModel <- FALSE

## Relative paths assuming the working directory is the script directory
## containing this script
scriptDir <- getwd()
projectDir <- scriptDir
figDir <- file.path(projectDir, "deliv", "figure")
tabDir <- file.path(projectDir, "deliv", "table")
dataDir <- file.path(projectDir, "data", "derived")
outDir <- file.path(projectDir, "output")
toolsDir <- file.path(scriptDir, "tools")
invisible(dir.create(figDir, recursive = TRUE))
invisible(dir.create(tabDir, recursive = TRUE))
invisible(dir.create(outDir, recursive = TRUE))

suppressMessages(library(bayesplot))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
library(cmdstanr)

set_cmdstan_path(path = file.path("..", "..", "cmdstan"))

set.seed(10271998)

nChains <- 4
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin

mod <- cmdstan_model(paste0(modelName, ".stan"),
                     quiet=FALSE,
                     compiler_flags = c("STANC2=true"))

## alternatively, for MPI applications
## mod <- cmdstan_model(paste0(modelName, ".stan"),
##                      quiet=FALSE,
##                      compiler_flags=c("TORSTEN_MPI=1",
##                                       "CXX=mpicxx",
##                                       "TBB_CXX_TYPE=clang",
##                                       "STANC2=true"))

initFiles <- paste0("init", 1:nChains, ".json")

fit <- mod$sample(data = "data.json", init=initFiles, seed = 11191951,
                  num_chains=nChains,
                  num_cores = min(nChains, detectCores()),
                  num_warmup = nBurnin,
                  num_samples = nIter - nBurnin,
                  thin = nThin,
                  refresh = 10,
                  output_dir = file.path(outDir),
                  adapt_delta = 0.95)
