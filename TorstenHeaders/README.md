# `TorstenHeaders` for `Torsten` library
This R package installs `Torsten` library. It installs the following dependencies

-  [rstan](https://github.com/stan-dev/rstan)
-  [stan for Torsten](https://github.com/stan-dev/stan)
-  [math for Torsten](https://github.com/stan-dev/math)

A few model examples can be found at Torsten's [example-models](https://github.com/metrumresearchgroup/example-models).

Installation
------------
```r
devtools::install_github('metrumresearchgroup/TorstenHeaders')
install_torsten()
```
or after cloning this repo
```r
devtools::load_all('TorstenHeaders')
install_torsten()
```
System requirement
------------------
See [Torsten](https://github.com/metrumresearchgroup/Torsten) for details.
