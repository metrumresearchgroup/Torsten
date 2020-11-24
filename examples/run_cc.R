## Aki's summary script
perf.cc <- function(stanfit) {
    (n_chain = stanfit@sim$chains)
    (n_warmup = stanfit@sim$warmup)
    n_iter = stanfit@sim$iter-n_warmup
    sampler_params <- rstan:::get_sampler_params(stanfit, inc_warmup = TRUE)
    leapfrogs = sapply(sampler_params, function(x) x[, "n_leapfrog__"])
    (sum_warmup_leapfrogs = sum(leapfrogs[1:n_warmup,]))
    (sum_leapfrogs = sum(leapfrogs[n_warmup+(1:n_iter),]))
    (mean_warmup_leapfrogs = sum_warmup_leapfrogs/n_warmup / n_chain)
    (mean_leapfrogs = sum_leapfrogs/n_iter / n_chain)
    mon = rstan::monitor(as.array(stanfit), warmup=0, print=FALSE)
    (maxrhat = max(mon[,'Rhat']))
    bulk_ess_per_iter = mon[,'Bulk_ESS']/n_iter / n_chain
    tail_ess_per_iter = mon[,'Tail_ESS']/n_iter / n_chain
    bulk_ess_per_leapfrog = mon[,'Bulk_ESS']/sum_leapfrogs
    tail_ess_per_leapfrog = mon[,'Tail_ESS']/sum_leapfrogs
    min(bulk_ess_per_iter)
    min(tail_ess_per_iter)
    min(bulk_ess_per_leapfrog)
    min(tail_ess_per_leapfrog)
    elapsed <- as.data.frame(rstan::get_elapsed_time(stanfit))
    (stepsizes = sapply(sampler_params, function(x) x[, "stepsize__"])[n_iter,])

    res <- data.frame(run = c(sum_warmup_leapfrogs / n_chain, sum_leapfrogs / n_chain,
                                   mean_warmup_leapfrogs, mean_leapfrogs,
                                   min(bulk_ess_per_iter),
                              min(tail_ess_per_iter),
                              min(bulk_ess_per_leapfrog),
                              min(tail_ess_per_leapfrog),
                              max(elapsed$warmup + elapsed$sample)))
    row.names(res) <- c("leapfrogs(warmup)", "leapfrogs(sampling)",
                        "leapfrogs(warmup)/iter", "leapfrogs(sampling)/iter",
                        "min(bulk_ESS/iter)", "min(tail_ESS/iter)",
                        "min(bulk_ESS/leapfrog)", "min(tail_ESS/leapfrog)",
                        "max(elapsed_time)")
    return(res)
}

## run mpi cross chain job
run.mpi <- function(modelpath, model, metric, np, hostfile, seed, init, adapt.arg="") {
    model.file = paste(model,"stan", sep="")
    data.file = paste(model,".data.R", sep="")
    system(paste("make -j4 ", file.path(modelpath, model, model), sep=""))

    cmdwd <- getwd();
    setwd(file.path(modelpath, model))

    rng.seed = sample(1:999999, 1)
    if (!missing(seed)) {
        rng.seed = seed
    }

    if (missing(init)) {
        ## system(paste("mpiexec -n 8 -l -bind-to core ./", model, " sample save_warmup=1 algorithm=hmc metric=", metric, " adapt ", adapt.arg, " data file=", data.file, " random seed=", seed, " ",sep=""))
        run.string <- paste("mpiexec -bind-to core -n ",np," -l -f ", hostfile, " ./", model, " sample save_warmup=1 algorithm=hmc metric=", metric, " adapt ", adapt.arg, " data file=", data.file, " random seed=", seed, " ",sep="")
        print(run.string)
        system(run.string)
        fit <- rstan::read_stan_csv(dir(pattern="mpi.[0-9]*.output.csv", full.name=TRUE))
        summary <- data.frame(perf.cc(fit))
    } else {
        run.string <- paste("mpiexec -bind-to core -n ",np," -l -f ", hostfile, " ./", model, " sample save_warmup=1 algorithm=hmc metric=", metric, " adapt ", adapt.arg, " data file=", data.file, " init=", init, " random seed=", seed, " ",sep="")
        print(run.string)
        system(run.string)
        fit <- rstan::read_stan_csv(dir(pattern="mpi.[0-9]*.output.csv", full.name=TRUE))
        summary <- data.frame(perf.cc(fit))
    }
    setwd(cmdwd)
    return(summary)
}

## here seed is a sequence, with each member fed to one MPI run
multiple.run.mpi <- function(modelpath, model, metric, np, hostfile, seed, init, adapt.arg="") {
    if (missing(init)) {
        res <- lapply(seed, function(s) {run.mpi(modelpath, model, metric, np, hostfile, s, adapt.arg=adapt.arg)})
    } else {
        res <- lapply(seed, function(s) {run.mpi(modelpath, model, metric, np, hostfile, s, init, adapt.arg=adapt.arg)})
    }
    summary <- do.call(cbind, res)
    n <- ncol(summary)
    summary$avg <- apply(summary[,1:n], 1, mean)
    summary$sd <- apply(summary[,1:n], 1, sd)
    return(summary)
}

## run mpi cross chain job
run.seq <- function(modelpath, model, metric, nchain, hostfile, seed, init) {
    model.file = paste(model,"stan", sep="")
    data.file = paste(model,".data.R", sep="")
    system(paste("make -j4 ", file.path(modelpath, model, model), sep=""))

    cmdwd <- getwd();
    setwd(file.path(modelpath, model))

    rng.seed = sample(1:999999, 1)
    if (!missing(seed)) {
        rng.seed = seed
    }

    mpi.string  <- paste("mpiexec -bind-to core -f ", hostfile)
    app.string.1 <- paste0(" -n 1 ./", model, "_seq", " sample save_warmup=1 algorithm=hmc metric=", metric, " data file=", data.file, " random seed=", seed)
    if (missing(init)) {
        app.string <- paste(sapply(1:nchain, function(chain) {
            paste0(app.string.1,
                   " output file=", paste0("seq.", chain, ".output.csv"),
                   " id=", chain) }), collapse=" :")
        print(paste0(mpi.string, app.string))
        system(paste0(mpi.string, app.string))
        fit <- rstan::read_stan_csv(dir(pattern="^seq.[0-9]*.output.csv", full.name=TRUE))
        summary <- data.frame(perf.cc(fit))
    } else {
        app.string <- paste(sapply(1:nchain, function(chain) {
            paste0(app.string.1, paste(" init=", "mpi.", chain, init),
                   " output file=", paste0("seq.", chain, ".output.csv"),
                   " id=", chain) }), collapse=" :")
        print(paste0(mpi.string, app.string))
        system(paste0(mpi.string, app.string))
        fit <- rstan::read_stan_csv(dir(pattern="^seq.[0-9]*.output.csv", full.name=TRUE))
        summary <- data.frame(perf.cc(fit))
    }
    setwd(cmdwd)
    return(summary)
}

## here seed is a sequence, with each member fed to one MPI run
multiple.run.seq <- function(modelpath, model, metric, np, hostfile, seed, init) {
    if (missing(init)) {
        res <- lapply(seed, function(s) {run.seq(modelpath, model, metric, np, hostfile, s)})
    } else {
        res <- lapply(seed, function(s) {run.seq(modelpath, model, metric, np, hostfile, s, init)})
    }
    summary <- do.call(cbind, res)
    n <- ncol(summary)
    summary$avg <- apply(summary[,1:n], 1, mean)
    summary$sd <- apply(summary[,1:n], 1, sd)
    return(summary)
}

multiple.run.summary <- function(modelpath, model, np, hostfile, seed, init, adapt.arg="") {
    res.mpi.diag <- multiple.run.mpi(modelpath, model, "diag_e", np, hostfile, seed, init, adapt.arg)
    res.mpi.diag$metric <- "diag_e"
    res.mpi.diag$Stan_run <- "cross_chain"
    res.seq.diag <- multiple.run.seq(modelpath, model, "diag_e", np, hostfile, seed, init)
    res.seq.diag$metric <- "diag_e"
    res.seq.diag$Stan_run <- "regular"

    res.mpi.dense <- multiple.run.mpi(modelpath, model, "dense_e", np, hostfile, seed, init, adapt.arg)
    res.mpi.dense$metric <- "dense_e"
    res.mpi.dense$Stan_run <- "cross_chain"
    res.seq.dense <- multiple.run.seq(modelpath, model, "dense_e", np, hostfile, seed, init)
    res.seq.dense$metric <- "dense_e"
    res.seq.dense$Stan_run <- "regular"

    res <- do.call(rbind, lapply(list(res.mpi.diag, res.seq.diag, res.mpi.dense, res.seq.dense),
                                 function(x) {y <- data.frame(row.names(x), x$avg, x$sd, x$metric, x$Stan_run);names(y) <- c("performance", "avg", "sd", "metric", "run");return(y)}))

    library("ggplot2")
    ggplot(res, aes(x=metric)) + geom_bar(aes(y=avg,fill=run),position=position_dodge(0.8),stat="identity",alpha=0.7,width=0.8) + geom_errorbar(aes(ymin=avg-sd,ymax=avg+sd,group=run),colour="black",alpha=0.4,size=0.4,width=0.3,position=position_dodge(0.8)) +
        facet_wrap(performance ~ ., scales="free_y")
    ggsave(file=file.path(modelpath, model, "cross_chain_summary.png"))

    return(res)
}

## effect of target ESS, "target.ess" is a sequence s.a c(100, 200, 400, 800)
multiple.run.ess <- function(modelpath, model, np, nchains, hostfile, seed, target.ess)
{
    n <- length(seed)
    res.mpi.diag <- lapply(target.ess,
                           function(ess) {
                               res  <- multiple.run.mpi(modelpath, model, "diag_e", np, hostfile, seed,
                                                        adapt.arg=paste0("num_cross_chains=",nchains," cross_chain_ess=",ess))
                               res$Stan_run <- paste0("cross chain: target ess=",ess)
                               res$metric <- "diag_e"
                               res <- data.frame(row.names(res), res$avg, res$sd, res$metric, res$Stan_run)
                               return(res)
                           })
    res.mpi.diag <- do.call(rbind,res.mpi.diag)
    names(res.mpi.diag) <- c("performance", "avg", "sd", "metric", "run")

    res.mpi.dense <- lapply(target.ess,
                           function(ess) {
                               res  <- multiple.run.mpi(modelpath, model, "dense_e", np, hostfile, seed,
                                                        adapt.arg=paste0("cross_chain_ess=",ess))
                               res$Stan_run <- paste0("cross chain: target ess=",ess)
                               res$metric <- "dense_e"
                               res <- data.frame(row.names(res), res$avg, res$sd, res$metric, res$Stan_run)
                               return(res)
                           })
    res.mpi.dense <- do.call(rbind,res.mpi.dense)
    names(res.mpi.dense) <- c("performance", "avg", "sd", "metric", "run")

    res.seq.diag <- multiple.run.seq(modelpath, model, "diag_e", np, hostfile, seed)
    res.seq.diag$metric <- "diag_e"
    res.seq.diag$Stan_run <- "regular"
    res.seq.diag <- data.frame(row.names(res.seq.diag), res.seq.diag$avg, res.seq.diag$sd, res.seq.diag$metric, res.seq.diag$Stan_run)
    names(res.seq.diag) <- c("performance", "avg", "sd", "metric", "run")

    res.seq.dense <- multiple.run.seq(modelpath, model, "dense_e", np, hostfile, seed)
    res.seq.dense$metric <- "dense_e"
    res.seq.dense$Stan_run <- "regular"
    res.seq.dense <- data.frame(row.names(res.seq.dense), res.seq.dense$avg, res.seq.dense$sd, res.seq.dense$metric, res.seq.dense$Stan_run)
    names(res.seq.dense) <- c("performance", "avg", "sd", "metric", "run")

    res <- rbind(res.mpi.diag, res.mpi.dense, res.seq.diag, res.seq.dense)

    library("ggplot2")
    ggplot(res, aes(x=metric)) + geom_bar(aes(y=avg,fill=run),position=position_dodge(0.8),stat="identity",alpha=0.7,width=0.8) + geom_errorbar(aes(ymin=avg-sd,ymax=avg+sd,group=run),colour="black",alpha=0.4,size=0.4,width=0.3,position=position_dodge(0.8)) +
        facet_wrap(performance ~ ., scales="free_y") +
        theme(legend.position="bottom") + guides(fill=guide_legend(nrow=2,byrow=TRUE))
    ggsave(file=file.path(modelpath, model, "cross_chain_ess_effect.png"))

    return(res)
}

## run both cross-chain & regular warmup version and output summary
run.cc.metric <- function(modelpath, model, metric, np, seed, init, adapt.arg="") {
    model.file = paste(model,"stan", sep="")
    data.file = paste(model,".data.R", sep="")
    system(paste("make -j4 ", file.path(modelpath, model, model), sep=""))
    cmdwd <- getwd();

    setwd(file.path(modelpath, model))

    if (missing(init)) {
        ## system(paste("mpiexec -n 8 -l -bind-to core ./", model, " sample save_warmup=1 algorithm=hmc metric=", metric, " adapt ", adapt.arg, " data file=", data.file, " random seed=", seed, " ",sep=""))
        system(paste("mpiexec -bind-to core -n",np,"-l ./", model, " sample save_warmup=1 algorithm=hmc metric=", metric, " adapt ", adapt.arg, " data file=", data.file, " random seed=", seed, " ",sep=""))
        fit.mpi <- rstan::read_stan_csv(dir(pattern="mpi.[0-3].output.csv", full.name=TRUE))
        system(paste("for i in {0..3}; do ./", model, " sample save_warmup=1 algorithm=hmc metric=", metric, " data file=", data.file, " random seed=", seed, " output file=seq.$i.output.csv id=$i;done", " ",sep=""))
        fit.seq <- rstan::read_stan_csv(dir(pattern="^seq.[0-3].output.csv", full.name=TRUE))
        summary <- data.frame(perf.cc(fit.mpi), perf.cc(fit.seq))
        colnames(summary) <- c(paste0("MPI.",metric), paste0("regular.",metric))        
    } else {
        system(paste("mpiexec -bind-to core -n",np,"-l ./", model, " sample save_warmup=1 algorithm=hmc metric=", metric, " adapt ", adapt.arg, " data file=", data.file, " init=", init, " random seed=", seed, " ",sep=""))
        fit.mpi <- rstan::read_stan_csv(dir(pattern="mpi.[0-3].output.csv", full.name=TRUE))
        system(paste("for i in {0..3}; do ./", model, " sample save_warmup=1 algorithm=hmc metric=", metric, " data file=", data.file, " init=mpi.$i.", init, " random seed=", seed, " output file=seq.$i.output.csv id=$i;done", " ",sep=""))
        fit.seq <- rstan::read_stan_csv(dir(pattern="^seq.[0-3].output.csv", full.name=TRUE))
        summary <- data.frame(perf.cc(fit.mpi), perf.cc(fit.seq))
        colnames(summary) <- c(paste0("MPI.",metric), paste0("regular.",metric))
    }

    setwd(cmdwd)
    return(list(mpi=fit.mpi, seq=fit.seq, summary=summary))
}

## run cross-chain & regular warmup for both diag & dense metric
run.cc <- function(modelpath, model, init, seed, adapt.arg="") {
    rng.seed = sample(1:999999, 1)
    if (!missing(seed)) {
        rng.seed = seed
    }

    ## seed = sample(999999:9999999, 1)

    fits.diag  <- run.cc.metric(modelpath, model, "diag_e",  rng.seed, init, adapt.arg)
    fits.dense <- run.cc.metric(modelpath, model, "dense_e", rng.seed, init, adapt.arg)

    library("gridExtra")
    library("bayesplot")
    library("dplyr")

    ## print out divergence diagno summary
    pdf(file.path(modelpath, model, "cross_chain_summary.pdf"), paper="a4")
    grid.table(format(cbind(fits.diag[["summary"]], fits.dense[["summary"]]), digits=3))
    lp <- log_posterior(fits.diag[["mpi"]])
    np <- nuts_params(fits.diag[["mpi"]])
    p <- grid.arrange(mcmc_nuts_divergence(np, lp),top="mcmc_nuts_divergence: cross-chain + diag_e")
    print(p)
    lp <- log_posterior(fits.diag[["seq"]])
    np <- nuts_params(fits.diag[["seq"]])
    p <- grid.arrange(mcmc_nuts_divergence(np, lp),top="mcmc_nuts_divergence: regular + diag_e")
    print(p)
    lp <- log_posterior(fits.dense[["mpi"]])
    np <- nuts_params(fits.dense[["mpi"]])
    p <- grid.arrange(mcmc_nuts_divergence(np, lp),top="mcmc_nuts_divergence: cross-chain + dense_e")
    print(p)
    lp <- log_posterior(fits.dense[["seq"]])
    np <- nuts_params(fits.dense[["seq"]])
    p <- grid.arrange(mcmc_nuts_divergence(np, lp),top="mcmc_nuts_divergence: regular + dense_e")
    print(p)
    dev.off()

    grid.table(format(cbind(fits.diag[["summary"]], fits.dense[["summary"]]), digits=3))

    d <- cbind(fits.diag[["summary"]], fits.dense[["summary"]])
    d%>%
        tibble::rownames_to_column()%>%
        tidyr::gather('col','val',-rowname)%>%
        tidyr::separate(col,c('type','metric'),sep='\\.')%>%
        dplyr::filter(grepl('',rowname))%>%
        ggplot2::ggplot(ggplot2::aes(x=metric,y=val,fill=type)) +
        ggplot2::geom_bar(stat = 'identity',position = ggplot2::position_dodge()) +
        ggplot2::facet_wrap(~rowname,scales='free',ncol=2)

    ggplot2::ggsave(file=file.path(modelpath, model, "cross_chain_summary.png"))

    return(list(summary=cbind(fits.diag[["summary"]], fits.dense[["summary"]]),
                mpi.diag=fits.diag[["mpi"]],
                seq.diag=fits.diag[["seq"]],
                mpi.dense=fits.dense[["mpi"]],
                seq.dense=fits.dense[["seq"]]))
}
