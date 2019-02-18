default_stan_parameters <- function() {         #parameters from #stan model
    par <- c("cHatDat1",
             "cHatDat2",
             "cHatDat3",
             "cHatDat4",
             "cHatDat5",
             "cHatDat6",
             "cHatDat7",
             "cHatDat8",
             "cHatObs1",
             "cHatObs2",
             "cHatObs3",
             "cHatObs4",
             "cHatObs5",
             "cHatObs6",
             "cHatObs7",
             "cHatObs8",
             "cHatDPP1",
             "cHatDPP2",
             "cHatDPP3",
             "cHatDPP4",
             "cHatDPP5",
             "cHatDPP6",
             "cHatDPP7",
             "cHatDPP8",
             "cHatPDP1",
             "cHatPDP2",
             "cHatPDP3",
             "cHatPDP4",
             "cHatPDP5",
             "cHatPDP6",
             "cHatPDP7",
             "cHatPDP8",
             "cHatPPD1",
             "cHatPPD2",
             "cHatPPD3",
             "cHatPPD4",
             "cHatPPD5",
             "cHatPPD6",
             "cHatPPD7",
             "cHatPPD8",
             "cHatDDP1",
             "cHatDDP2",
             "cHatDDP3",
             "cHatDDP4",
             "cHatDDP5",
             "cHatDDP6",
             "cHatDDP7",
             "cHatDDP8",
             "cHatDPD1",
             "cHatDPD2",
             "cHatDPD3",
             "cHatDPD4",
             "cHatDPD5",
             "cHatDPD6",
             "cHatDPD7",
             "cHatDPD8",
             "cHatPDD1",
             "cHatPDD2",
             "cHatPDD3",
             "cHatPDD4",
             "cHatPDD5",
             "cHatPDD6",
             "cHatPDD7",
             "cHatPDD8"
             )
    return (par)
}

fit.from.data <- function(data,
                          modelName,
                          stan_parameters,
                          endtime = tail(data$time, n=1) + 1.0) {

torstenDir <- Sys.getenv("TORSTEN_PATH")
try(if(torstenDir == "") stop("Missing environment variable 'TORSTEN_PATH'"))
modelDir <- file.path(torstenDir, "tests")

library(rstan)
library(parallel)
library(mrgsolve)  ## tested with version 0.7.6.9029
library(testthat)

set.seed(11091989) ## not required but assures repeatable results

    rstan_options(auto_write = TRUE)

    parametersToPlot <- stan_parameters
    parameters <- c(parametersToPlot)
    parametersToPlot <- c("lp__", parametersToPlot)

    nChains <- 1 ## don't really need this
    nIter <- 1;
    fit <- stan(file = file.path(modelDir, paste0(modelName, ".stan")),
                data = data,
                pars = parameters,
                iter = nIter,
                chains = nChains,
                algorithm="Fixed_param",
                cores = min(nChains, parallel::detectCores()))
    return(fit)
}

simulation_test <- function(event,
                            steptime,
                            modelName,
                            mrgcode,
                            stan_parameters,
                            need.fit = TRUE,
                            endtime = tail(steptime, n=1) + 1.0) {

torstenDir <- Sys.getenv("TORSTEN_PATH")
try(if(torstenDir == "") stop("Missing environment variable 'TORSTEN_PATH'"))
modelDir <- file.path(torstenDir, "tests")

library(rstan)
library(parallel)
library(mrgsolve)  ## tested with version 0.7.6.9029
library(testthat)

set.seed(11091989) ## not required but assures repeatable results

###############################################################################################
## Simulate deterministic data using mrgsolve

code <- mrgcode
mod <- mcode("accum", code) %>% Req(DV) %>% update(end=endtime,delta=0.1)
e1 <- event
mod %>% ev(e1) %>% mrgsim(end = endtime) %>% plot # plot data
time <- steptime
SimData <-
  mod %>%
  ev(e1) %>%
  carry_out(cmt, ii, addl, rate, amt, evid, ss) %>%
  mrgsim(Req = "DV", end = -1, add = time, recsort = 3) %>%
  as.data.frame
## print(SimData)
head(SimData)
SimData$cmt[SimData$cmt == 0] <- 2 ## adjust cmt (adopt NONMEM convention)
SimData <- SimData[!((SimData$evid == 0)&(SimData$DV == 0)),] ## remove observation with 0 drug concentration

################################################################################################
# Format Data for Stan using RStan
nt <- nrow(SimData)
iObs <- with(SimData, (1:nrow(SimData))[evid == 0])
nObs <- length(iObs)
data <- with(SimData,  ## create Stan data set
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

    if (need.fit) {               #stan fit & #mrgsolve data
        rstan_options(auto_write = TRUE)

        parametersToPlot <- stan_parameters
        parameters <- c(parametersToPlot)
        parametersToPlot <- c("lp__", parametersToPlot)

        nChains <- 1 ## don't really need this
        nIter <- 1;
        fit <- stan(file = file.path(modelDir, paste0(modelName, ".stan")),
                    data = data,
                    pars = parameters,
                    iter = nIter,
                    chains = nChains,
                    algorithm="Fixed_param",
                    cores = min(nChains, parallel::detectCores()))
        return(list("fit"=fit, "data"=data))
    } else {
        return (data)                    #only mrgsolve data
    }
}

do_tests <- function(fit, data, pars=default_stan_parameters(), tol=1e-6) {
    lapply(pars, (function(par) {
        ## ptable <- parameterTable(fit, par)
        ## only compare the 1st and last result in time
        ## expect_equal(data$cObs[1], ptable[1, "mean"]                    , tolerance=tol)
        ## expect_equal(tail(data$cObs, n=1), tail(ptable, n=1)[1, "mean"] , tolerance=tol)
        fit.res <- as.array(extract(fit, pars=par)[[par]])
        n <- length(fit.res)            #nb. of time steps
        for (i in 1:length(fit.res)) {
            expect_equal(data$cObs[i], fit.res[i], tolerance=tol)
        }
    }
    ))
}
