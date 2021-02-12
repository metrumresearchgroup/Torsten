modelName <- "pk3cpt_linode"
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
toolsDir <- file.path(dirname(projectDir), "R", "tools")
invisible(dir.create(figDir, recursive = TRUE))
invisible(dir.create(tabDir, recursive = TRUE))
invisible(dir.create(outDir, recursive = TRUE))

suppressMessages(library(tidyverse))
suppressMessages(library(bayesplot))
suppressMessages(library(parallel))
library(cmdstanr)

set_cmdstan_path(path = file.path(dirname(dirname(scriptDir)), "cmdstan"))

## reproducibility
set.seed(10271998)

nChains <- 4
nPost <- 500 ## Number of post-burn-in samples per chain after thinning
nBurn <- 500 ## Number of burn-in samples per chain after thinning
nThin <- 1

nIter <- (nPost + nBurn) * nThin
nBurnin <- nBurn * nThin
nsample = nIter - nBurnin

## compile model
mod <- cmdstan_model(paste0(modelName, ".stan"), quiet=FALSE)

## run MCMC
fit <- mod$sample(data = "pk3cpt_linode.data.R", init="pk3cpt_linode.init.R", seed = 3191951,
                  chains=nChains,
                  parallel_chains = min(nChains, detectCores()),
                  iter_warmup = nBurnin,
                  iter_sampling = nsample,
                  thin = nThin,
                  output_dir = file.path(outDir))

## figures
fit.summary <- fit$summary()
fit.summary %>%
    select(variable,rhat) %>%
    slice(1:7) %>%
    pivot_wider(names_from = variable,values_from = rhat) %>%
    unlist() %>%
    mcmc_rhat() + yaxis_text(hjust = 0)
ggsave(file.path(figDir, "diag_rhat.pdf"))

fit.summary %>%
    select(variable,ess_bulk) %>%
    slice(1:7) %>%
    pivot_wider(names_from = variable,values_from = ess_bulk) %>%
    `/`(nsample * nChains) %>%
    unlist() %>%
    mcmc_neff() + yaxis_text(hjust = 0)
ggsave(file.path(figDir, "diag_neff_bulk_ratio.pdf"))

fit.summary %>%
    select(variable,ess_tail) %>%
    slice(1:7) %>%
    pivot_wider(names_from = variable,values_from = ess_tail) %>%
    `/`(nsample * nChains) %>%
    unlist() %>%
    mcmc_neff() + yaxis_text(hjust = 0)
ggsave(file.path(figDir, "diag_neff_tail_ratio.pdf"))

p <- mcmc_trace(fit$draws(),  pars = c("CL", "Q", "V1", "V2", "ka", "sigma"), n_warmup = 0,
                facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 12)
ggsave(file.path(figDir, "history.pdf"))

mcmc_dens_overlay(fit$draws(), pars = c("CL", "Q", "V1", "V2", "ka", "sigma"),
                  facet_args = list(nrow = 2)) +
                  facet_text(size = 14)
ggsave(file.path(figDir, "density.pdf"))

p <- mcmc_pairs(fit$draws(), pars = c("CL", "Q", "V1", "V2", "ka", "sigma"),
                off_diag_args = list(size = 1, alpha = 0.5))
ggsave(file.path(figDir, "pair.pdf"), p, width = 10, height = 8)
