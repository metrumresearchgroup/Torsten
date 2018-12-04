context("GenTwoCptModel_rk45")

library(mrgsolve)  ## tested with version 0.7.6.9029
library(testthat)

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

run_test_GenTwoCptModel_rk45 <- function(event, steptime) {
modelName <- "test_GenTwoCptModel_rk45"
par <- default_stan_parameters()
return (simulation_test (event, steptime, modelName, mrgcode, par))
}

get_mrgdata_GenTwoCptModel_rk45 <- function(event, steptime) {
modelName <- "test_GenTwoCptModel_rk45"
par <- default_stan_parameters()
return (simulation_test (event, steptime, modelName, mrgcode, par, need.fit=FALSE))
}

testmodel = "test_GenTwoCptModel_rk45"
dose = "instantaneous bolus doses"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000)
    time <- c(0, 0.083, 0.167) ## Observation times
    res <- run_test_GenTwoCptModel_rk45(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000, rate = 8000)
    time <- c(0, 0.083, 0.167) ## Observation times
    res <- run_test_GenTwoCptModel_rk45(e1, time)
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
    res <- run_test_GenTwoCptModel_rk45(ev(c(e1, e2)), time)
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
    res <- run_test_GenTwoCptModel_rk45(ev(c(e1, e2)), time)
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
    res <- run_test_GenTwoCptModel_rk45(ev(c(e1, e2)), time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "steady state with multiple bolus"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 80 * 1000, rate = 0, ss=1, ii=12, addl=2)
    time <- c(0, 12, 22) ## Observation times
    res <- run_test_GenTwoCptModel_rk45(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "steady state with multiple infusion"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 10, rate=5, ss=1, ii=12, addl=2)
    time <- c(0, 12, 22) ## Observation times
    res <- run_test_GenTwoCptModel_rk45(e1, time)
    fit <- res$fit
    data <- res$data
    do_tests(fit, data)
})

dose = "steady state: bolus @ GUT and infusion @ CENT"
context(paste(testmodel, "with", dose))
test_that(dose, {
    e1 <- ev(amt = 800, rate = 800, ss=1, ii=12, addl=5, cmt=2)
    e2 <- ev(time=2, amt = 800, ii=6, addl=5)
    time <- c(0, 0.5, 1.0, 2.0, 3.0, 50.0) ## Observation times
    res <- run_test_GenTwoCptModel_rk45(e1+e2, time)
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
    modelName <- "test_GenTwoCptModel_rk45"
    par <- default_stan_parameters()
    data <- get_mrgdata_GenTwoCptModel_rk45(e1, time)
    data$amt[1] = 0
    data$ii[1] = 0
    data$rate = rep(800, data$nt)
    data$ss = rep(1, data$nt)
    data$addl = rep(0, data$nt)
    data$evid = rep(1, data$nt)
    fit <- fit.from.data(data, "test_GenTwoCptModel_rk45", par)
    do_tests(fit, data)
})

