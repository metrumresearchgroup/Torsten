
#' @title Installation for torsten
#' @description installation of torsten with rstan, that replaces stanHeaders.
#' @param StanHeaders_version package_version, package version of StanHeaders to append Torsten Default: NULL
#' @param rstan_version package_version, package version of rstan to install, Default: NULL
#' @param lib character, giving the library directory where to install the packages, Default: .libPaths()[1]
#' @param ... parameters to pass to install.packages
#' @details installation will replace the 'StanHeaders/include/src/stan' and 'StanHeaders/include/stan' of stanHeaders and install
#' source rstan without dependencies.
#' @rdname install_torsten
#' @importFrom devtools install_version
#' @importFrom utils install.packages installed.packages packageVersion
#' @export

install_torsten <- function(StanHeaders_version=NULL,
                            rstan_version=NULL,
                            lib=.libPaths()[1],
                            ...) {

  lib <- normalizePath(lib,winslash = '/')

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


  TH <- find.package('torstenHeaders')

  file.copy(file.path(TH,'stan'),file.path(lib,'StanHeaders/include/src/'),overwrite=TRUE,recursive=TRUE)
  file.copy(file.path(TH,'math/stan'),file.path(lib,'StanHeaders/include/'),overwrite=TRUE,recursive=TRUE)

  if(is.null(rstan_version)) rstan_version <- read.dcf(system.file('CURRENT_VERSION',package = 'torstenHeaders'),fields = 'rstan')

  rstan_version <- as.package_version(rstan_version)

  ## if(pkgVersionCRAN('rstan')==rstan_version){
  ##   utils::install.packages('rstan', lib=lib, type='source', ...)
  ## }else{
  ##   devtools::install_version(package = 'rstan', version = rstan_version, lib=lib, type='source', ...)
  ## }
  remotes::install_github('metrumresearchgroup/rstan/rstan/rstan', ref="torsten-develop")


}
