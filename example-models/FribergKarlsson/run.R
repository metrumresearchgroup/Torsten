library("cmdstanr")
library("bayesplot")
library("posterior")
library("ggplot2")
library("tidyr")
library("dplyr")

## use mpiexec for within-parapllel runs.
## mpiexec -n 8 -hostfile hostfile_ac -bind-to core ./FribergKarlsson sample num_warmup=700 num_samples=1000 data file=short_data.json init=init.R random seed=8765 id=3 output refresh=5 file=output.3.csv

## For within-chain parallel run, we do it in cmdstan. Here we only
## import the output.
fit <- as_cmdstan_fit(dir(pattern="output.[1-4].csv",full.names=TRUE))

pars <- c("CLHat"   ,  "QHat"     , "V1Hat"   ,  "V2Hat"   ,  "kaHat"   ,  "mttHat"   ,
"circ0Hat" ,  "alphaHat", "omega[1]" , "omega[2]",  "omega[3]",  "omega[4]", "omega[5]",
"omega[6]" ,  "omega[7]", "gamma",
"sigma"   ,  "sigmaNeut")

subset.pars <- subset_draws(fit$draws(), variable=pars)

## Write summary
write.csv(summarise_draws(subset.pars), file="summary_pars.csv")

## density plot
mcmc_dens_overlay(subset.pars, facet_args=list(ncol=4))+facet_text(size=10)+theme(axis.text=element_text(size=10))

ggsave("density.pdf", width=8, height=6)

## get data
data <- jsonlite::read_json("short_data.json")

## PPC PK
cobs.rep <-
as_draws_df(fit$draws(variables=c("cHatObsPred"))) %>%
    select(starts_with("cHatObsPred")) %>% as.matrix()

time <- data[["time"]]
nSubjects <- data[["nSubjects"]]
nObsPK  <- data[["nObsPK"]]
nObsPD  <- data[["nObsPD"]]
iObsPK  <- data[["iObsPK"]]
iObsPD  <- data[["iObsPD"]]

ppc_ribbon_grouped(y=data[["cObs"]], yrep=cobs.rep, x=time[iObsPK],
                   group=as.vector(sapply(1:nSubjects, function(i){rep(i, nObsPK/nSubjects)}))) + scale_x_continuous(name="time (h)") +
scale_y_continuous(name="drug plasma concentration (ng/mL)") + theme(axis.text=element_text(size=10))

ggsave("ppc_pk.pdf", width=8, height=6)

## PPC PD
neut.rep <- 
as_draws_df(fit$draws(variables=c("neutHatObsPred"))) %>%
    select(starts_with("neutHatObsPred")) %>% as.matrix()

ppc_ribbon_grouped(y=data[["neutObs"]], yrep=neut.rep, x=time[iObsPD],
                   group=as.vector(sapply(1:nSubjects, function(i){rep(i, nObsPD/nSubjects)}))) + scale_x_continuous(name="time (h)") +
scale_y_continuous(name="ANC") + theme(axis.text=element_text(size=10))

ggsave("ppc_pd.pdf", width=8, height=6)
