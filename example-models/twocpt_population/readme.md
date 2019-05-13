- [Description](#orgbd288cf)
- [Build](#org22c3eac)
  - [Edit/Add `cmdstan/make/local`](#org29043ee)
  - [Build in `cmdstan`](#org72f0a07)
- [Run](#org98b35be)
- [Results](#org0848f3a)


<a id="orgbd288cf"></a>

# Description

Parameter inference for the two-compartment solved using numerical integrator `pmx_integrate_ode_group_bdf`.


<a id="org22c3eac"></a>

# Build


<a id="org29043ee"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org72f0a07"></a>

## Build in `cmdstan`

```sh
make ../example-models/twocpt_population/twocpt_population
```


<a id="org98b35be"></a>

# Run

```sh
mpiexec -n 2 twocpt_population sample data file=twocpt_population.data.R init=twocpt_population.init.R
```


<a id="org0848f3a"></a>

# Results

MPI run results with 1, 2, 4, 8 processes, respectively, for a population size 8.

| run                | wall time(s) |
|------------------ |------------ |
| n<sub>proc</sub>=1 | 10274        |
| n<sub>proc</sub>=2 | 5411         |
| n<sub>proc</sub>=4 | 3043         |
| n<sub>proc</sub>=8 | 2654         |
