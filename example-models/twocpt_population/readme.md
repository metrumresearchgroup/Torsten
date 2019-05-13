- [Description](#orgc398e7c)
- [Build](#org4656f45)
  - [Edit/Add `cmdstan/make/local`](#orgc033e9e)
  - [Build in `cmdstan`](#org80e506d)
- [Run](#orgcc9ca1c)
- [Results](#org1bb727e)


<a id="orgc398e7c"></a>

# Description

Parameter inference for the two-compartment solved using numerical integrator `pmx_integrate_ode_group_bdf`.


<a id="org4656f45"></a>

# Build


<a id="orgc033e9e"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
CXXFLAGS += -isystem /usr/local/include                 # path to MPI library's headers
```


<a id="org80e506d"></a>

## Build in `cmdstan`

```sh
make ../example-models/twocpt_population/twocpt_population
```


<a id="orgcc9ca1c"></a>

# Run

```sh
mpiexec -n 2 twocpt_population sample data file=twocpt_population.data.R init=twocpt_population.init.R
```


<a id="org1bb727e"></a>

# Results

MPI run results with 1, 2, 4, 8 processes, respectively, for a population size 8.

| run                | wall time(s) |
|------------------ |------------ |
| n<sub>proc</sub>=1 | 10274        |
| n<sub>proc</sub>=2 | 5411         |
| n<sub>proc</sub>=4 | 3043         |
| n<sub>proc</sub>=8 | 2654         |
