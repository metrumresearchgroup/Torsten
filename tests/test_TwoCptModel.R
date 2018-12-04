context("TwoCptModel")

library(mrgsolve)  ## tested with version 0.7.6.9029
library(testthat)
testmodel = "TwoCptModel"

torstenDir <- Sys.getenv("TORSTEN_PATH")
try(if(torstenDir == "") stop("Missing environment variable 'TORSTEN_PATH'"))
source(file.path(torstenDir, "tests", "simulation_test_template.R"))

mrgcode <- '
$PARAM CL = 5, Q = 8, V2 = 20, V3 = 70, KA = 1.2

$PKMODEL cmt="GUT CENT PERI", depot = TRUE

$GLOBAL
#define CP (CENT/V2)

$TABLE
capture DV = CP * exp(EPS(1));

$CAPTURE CP
'

run_test_TwoCptModel <- function(event, steptime) {
modelName <- "test_TwoCptModel"
par <- default_stan_parameters()
return (simulation_test (event, steptime, modelName, mrgcode, par))
}

get_mrgdata_TwoCptModel <- function(event, steptime) {
modelName <- "test_TwoCptModel"
par <- default_stan_parameters()
return (simulation_test (event, steptime, modelName, mrgcode, par, need.fit=FALSE))
}

dose = "instantaneous bolus doses"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000)
    time <- c(0, 0.083, 0.167) ## Observation times
    res <- run_test_TwoCptModel(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 8000, rate = 8000)
    time <- c(0, 0.083, 0.167, 2.0) ## Observation times
    res <- run_test_TwoCptModel(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "overlapping infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(time=0.0, amt = 80 * 1000, rate = 8000)
    e2 <- ev(time=6.0, amt = 80 * 1000, rate = 4000)
    time <- c(0, 5, 10) ## Observation times
    res <- run_test_TwoCptModel(ev(c(e1, e2)), time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "bolus + infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(time=0.0, amt = 80 * 1000)
    e2 <- ev(time=1.0, amt = 80 * 1000, rate = 8000)
    time <- c(0, 0.083, 0.167, 6.0) ## Observation times
    res <- run_test_TwoCptModel(ev(c(e1, e2)), time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "multiple bolus + infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(time=0.0, amt = 80 * 1000, ii = 6, addl = 4)
    e2 <- ev(time=1.0, amt = 80000, rate = 10000, ii = 6, addl=4)
    time <- c(0, 18.0, 36.0) ## Observation times
    res <- run_test_TwoCptModel(ev(c(e1, e2)), time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "steady state with multiple bolus"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000, rate = 0, ss=1, ii=12, addl=2)
    time <- c(0, 15.0, 22.0) ## Observation times
    res <- run_test_TwoCptModel(e1, time)
    fit <- res$fit
    data <- res$data
    lapply(default_stan_parameters(), (function(par) {
        fit.res <- as.array(extract(fit, pars=par)[[par]])
        n <- length(fit.res)            #nb. of time steps
        expect_equal(data$cObs[1], fit.res[1]                  , tolerance=1.e-3, scale=data$cObs[1])
        expect_equal(data$cObs[n], fit.res[n]                  , tolerance=1.e-3, scale=data$cObs[2])

        ## ptable <- parameterTable(fit, par)
        ## expect_equal(data$cObs[1], ptable[1, "mean"], tolerance=1.e-3, scale=data$cObs[1])
        ## expect_equal(data$cObs[2], ptable[2, "mean"], tolerance=1.e-3, scale=data$cObs[2])
        }
    ))
})

dose = "steady state with multiple infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 10, rate=5, ss=1, ii=12, addl=2)
    time <- c(0, 12, 22) ## Observation times
    res <- run_test_TwoCptModel(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "steady state with dosing at two compartments"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 800, rate = 800, ss=1, ii=6, addl=5, cmt=2)
    e2 <- ev(time=0.5, amt = 800, ii=6, addl=5)
    time <- c(0, 0.5, 1.0, 2.0, 3.0, 24.0) ## Observation times
    res <- run_test_TwoCptModel(e1+e2, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

## mrgsolve cannot do NONMEN's steady state with constant
## infusion, so we use a workaground here. The Torsten input
## is NONMEN style but the mrgsolve event is not.
dose = "steady state with const infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 700, rate = 800, ss=1, ii=700/800, addl=5, cmt=2)
    time <- c(0.0, 4.0) ## Observation times
    modelName <- "test_TwoCptModel"
    par <- default_stan_parameters()
    data <- get_mrgdata_TwoCptModel(e1, time)
    data$amt[1] = 0
    data$ii[1] = 0
    data$rate = rep(800, data$nt)
    data$ss = rep(1, data$nt)
    data$addl = rep(0, data$nt)
    data$evid = rep(1, data$nt)
    fit <- fit.from.data(data, "test_TwoCptModel", par)
    do_tests(fit, data)
})

