library("rstan")
library("parallel")

set.seed(11091989) ## not required but assures repeatable results

rstan_options(auto_write = TRUE)

n <- 2L
N_subj <- 3L
nt <- 7L
ts <- c(0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 100.0,
        0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 100.0,
        0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 100.0)
y1_obs <- c(0.165673,-0.242994,0.916383,1.04467,1.89353,-1.4577,1.24554,
            -0.194925,-1.20269,-0.974754,-0.792291,-0.479865,-1.0271,1.37066,
            -0.194925,-1.20269,-0.974754,-0.792291,-0.479865,-1.0271,1.37066)
y2_obs <- c(1.57484,-0.143474,0.593943,-0.733714,-1.24777,-1.50113,-0.0907854,
            0.464024,0.0336547,0.181032,0.0712328,1.12699,1.06031,0.338586,
            0.464024,0.0336547,0.181032,0.0712328,1.12699,1.06031,0.338586)

data <- list(n=n, N_subj=N_subj, nt=nt, ts=ts, y1_obs=y1_obs, y2_obs=y2_obs)

Sys.setenv(PKG_CXXFLAGS = "-DTORSTEN_MPI")

nChains <- 1 ## don't really need this
nIter <- 1000;
fit <- stan(file = "sho_group.stan",
            data = data,
            iter = nIter,
            chains = nChains,
            algorithm="NUTS",
            ## verbose=TRUE,
            cores = 1)

