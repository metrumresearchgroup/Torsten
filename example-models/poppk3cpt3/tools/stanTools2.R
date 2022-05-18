mcmcHistory <- function(posterior, pars = dimnames(posterior)[[3]], nParPerPage = 6, myTheme = NULL){
  ## Create MCMC history plots
  ## posterior = 3-D array of MCMC samples. Dims = {iterations, chains, parameters}. 
  ## pars = names of parameters to plot
  ## nParPerPage = maximum number of parameters to plot per page
  ## myTheme = ggplot2 theme
  
  require(bayesplot)
  posterior <- posterior[, , unlist(sapply(c(paste("^", pars, "$",sep=""),
                                             paste("^", pars, "\\[",sep="")),
                                           grep, 
                                           x = dimnames(posterior)[[3]]))]
  pars <- dimnames(posterior)[[3]]
  nPars <- length(pars)
  nPages <- ceiling(nPars / nParPerPage)
  parameters <- data.frame(parameter = pars,
                           page = sort(rep(1:nPages, length = nPars)),
                           stringsAsFactors = FALSE)
  
  lapply(1:nPages,
         function(i){
           postPage <- posterior[, , with(parameters, pars[page == i])]
           res <- mcmc_trace(postPage,
                             facet_args = list(ncol = 1, strip.position = "left")) +
             myTheme +
             scale_x_continuous(breaks = seq(0, nPost, len = 5))
           ## Code to remedy an apparent bug in Bayesplot
           res$data$Value <- res$data$value
           res
         }
  )
}

mcmcDensity <- function(posterior, pars = dimnames(posterior)[[3]], byChain = FALSE, 
                        nParPerPage = 16, myTheme = NULL, prior = NULL){
  ## Create density plots for marginal distributions of MCMC samples
  ## posterior = 3-D array of MCMC samples. Dims = {iterations, chains, parameters}
  ## pars = names of parameters to plot
  ## byChain = logical indicating whether to plot sensities by chain
  ## nParPerPage = maximum number of parameters to plot per page
  ## myTheme = ggplot2 theme
  ## prior = data.frame with columns value, density and 
  require(bayesplot)
  posterior <- posterior[, , unlist(sapply(c(paste("^", pars, "$",sep=""),
                                             paste("^", pars, "\\[",sep="")),
                                           grep, 
                                           x = dimnames(posterior)[[3]]))]
  pars <- dimnames(posterior)[[3]]
  
  nPars <- length(pars)
  nPages <- ceiling(nPars / nParPerPage)
  parameters <- data.frame(parameter = pars,
                           page = sort(rep(1:nPages, length = nPars)),
                           stringsAsFactors = FALSE)
  
  if(!is.null(prior)) prior <- prior %>% left_join(parameters)
  
  lapply(1:nPages,
         function(i){
           postPage <- posterior[, , with(parameters, pars[page == i])]
           if(byChain){
             p1 <- mcmc_dens_overlay(postPage)
           }else{
             p1 <- mcmc_dens(postPage)
           }
           if(!is.null(prior))
             p1 <- p1 + geom_line(data = subset(prior, page == i), 
                                  aes(x = value, y = density),
                                  color = "red")
           p1 <- p1 + myTheme
           ## Code to remedy an apparent bug in Bayesplot
           #p1$data$Value <- p1$data$value
           p1
         })
}


parameterTable <- function(posterior, pars = dimnames(posterior)[[3]]){
  require(rstan)
  posterior <- posterior[, , unlist(sapply(c(paste("^", pars, "$",sep=""),
                                             paste("^", pars, "\\[",sep="")),
                                           grep, 
                                           x = dimnames(posterior)[[3]]))]
  rstan:::monitor(posterior, warmup = 0, print = FALSE)
}
