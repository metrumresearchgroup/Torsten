## Install RStan, with Torsten built inside of it.
## WARNING: This file won't installs RStan's dependencies.

## Adjust directories to your setting.
scriptDir <- getwd()
lib <- file.path(scriptDir, "lib")
.libPaths(lib)

## Download and edit StanHeaders (version 2.16)
install.packages("https://cran.r-project.org/src/contrib/StanHeaders_2.16.0-1.tar.gz", repos = NULL, type = "source")
setwd(lib)
system("git clone https://github.com/metrumresearchgroup/stan.git")
setwd("stan")
# system("git checkout torsten-develop")  # comment out to get torsten-master
setwd(lib)
system("rm -rf StanHeaders/include/src/stan")
system("mv stan/src/stan StanHeaders/include/src/stan")
system("rm -rf stan")
system("git clone https://github.com/metrumresearchgroup/math.git")
setwd("math")
# system("git checkout torsten-develop")  # comment out to get torsten-master
setwd(lib)
system("rm -rf StanHeaders/include/stan")
system("mv math/stan StanHeaders/include/stan")
system("rm -rf math")

## Download rstan 2.16.2 without the dependencies.
install.packages("https://cran.r-project.org/src/contrib/rstan_2.16.2.tar.gz", repos = NULL, type = "source")

setwd(scriptDir)
