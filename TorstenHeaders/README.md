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
Please ensure the R toolchain includes a C++ compiler with C++11 support. In particular, R 3.3.0 and later is recommended as it contains toolchain based on gcc 4.9.3. On Windows platform, this further implies Rtools33 and later is recommended. 

Please ensure `CXXFLAGS` in `.R/Makevars` constains flag `-std=c++11`.
