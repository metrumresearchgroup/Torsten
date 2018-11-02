#' @title Installation for torsten through remote repository
#' @description installation of torsten with rstan, that replaces stanHeaders.
#' This function allows to specify which branch of torsten to install from.
#' @param StanHeaders_version package_version, package version of StanHeaders to append Torsten Default: NULL
#' @param rstan_version package_version, package version of rstan to install, Default: NULL
#' @param math_branch character, install the current build ('master') or different
#' development branch of torsten math library (e.g. 'develop'), Default: 'torsten-master'
#' @param stan_branch character, install the current build ('master') or different
#' development branch of torsten stan library (e.g. 'develop'), Default: 'torsten-master'
#' @param lib character, giving the library directory where to install the packages, Default: .libPaths()[1]
#' @param ... parameters to pass to install.packages
#' @details installation will replace the 'StanHeaders/include/src/stan' and 'StanHeaders/include/stan' of stanHeaders and install
#' source rstan without dependencies.
#' @rdname install_torsten_remote
#' @importFrom devtools install_version
#' @importFrom utils install.packages installed.packages packageVersion
#' @export

install_torsten_remote <- function(
  StanHeaders_version=NULL,
  rstan_version=NULL,
  math_branch='torsten-master',
  stan_branch='torsten-master',
  lib=.libPaths()[1],
  ...) {

  lib <- normalizePath(lib,winslash = '/')

  thiswd <- getwd()

  install_headers <- FALSE

  if(is.null(StanHeaders_version)){
    if(c('StanHeaders')%in%row.names(utils::installed.packages(lib.loc = lib))){
      StanHeaders_version <- utils::packageVersion('StanHeaders',lib.loc = lib)
    }else{
      StanHeaders_version <- read.dcf(system.file('CURRENT_VERSION',package = 'torstenHeaders'),fields = 'StanHeaders')
      StanHeaders_version <- as.package_version(StanHeaders_version)
      install_headers <- TRUE
    }
  }

  if(install_headers){
    if(pkgVersionCRAN('StanHeaders')==StanHeaders_version){
      utils::install.packages('StanHeaders', lib=lib, ...)
    }else{
      devtools::install_version(package = 'StanHeaders',version = StanHeaders_version,lib=lib, ...)
    }
  }

  td <- file.path(tempdir(),'torsten')
  dir.create(td,showWarnings = FALSE)

  for(repo in c('stan','math')){

    branch <- eval(parse(text = sprintf('%s_branch',repo)))

    system(sprintf("git clone --depth 1 -b %s https://github.com/metrumresearchgroup/%s.git %s/%s",branch, repo,td,repo))

    }

    system(sprintf("rm -rf %s/StanHeaders/include/src/stan",lib))
    system(sprintf("mv %s/stan/src/stan %s/StanHeaders/include/src/",td,lib))

    system(sprintf("rm -rf %s/StanHeaders/include/stan",lib))
    system(sprintf("mv %s/math/stan %s/StanHeaders/include/",td,lib))

  if(is.null(rstan_version)) rstan_version <- read.dcf(system.file('CURRENT_VERSION',package = 'torstenHeaders'),fields = 'rstan')

  rstan_version <- as.package_version(rstan_version)

  if(pkgVersionCRAN('rstan')==rstan_version){
    utils::install.packages('rstan', lib=lib, type='source', ...)
  }else{
    devtools::install_version(package = 'rstan', version = rstan_version, lib=lib, type='source', ...)
  }

  unlink(file.path(td,'math'),recursive = TRUE,force=TRUE)
  unlink(file.path(td,'stan'),recursive = TRUE,force=TRUE)

}
