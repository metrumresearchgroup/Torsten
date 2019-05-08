- [Description](#org74e8041)
- [Build](#org40bba43)
  - [Edit/Add `cmdstan/make/local`](#org1293759)
  - [Build in `cmdstan`](#org4cd239b)
- [Run](#org6eebecc)
- [Results](#org605b785)


<a id="org74e8041"></a>

# Description

Parameter inference for the harmonic oscillator ODE using Torsten's ODE group integrator `pmx_integrate_ode_group_rk45`.


<a id="org40bba43"></a>

# Build


<a id="org1293759"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org4cd239b"></a>

## Build in `cmdstan`

```sh
make ../example-models/harmonic_oscillator_ode_group_model/sho_group
```


<a id="org6eebecc"></a>

# Run

```sh
mpiexec -n 2 sho_group sample data file=sho_group.data.R
```


<a id="org605b785"></a>

# Results

Three chains are run using

-   sequential run using Stan's `bdf` integrator(output `sample1-3.stan.csv`)
-   sequential run using Torsten's `bdf` integrator(output `sample1-3.pmx.csv`)
-   MPI run using Torsten's `bdf` group integrator(output `sample1-3.mpi.csv`)

| chain | sequential stan | sequential pmx | MPI pmx |
|----- |--------------- |-------------- |------- |
| 1     | 31              | 21             | 10      |
| 2     | 31              | 21             | 10      |
| 3     | 32              | 21             | 10      |
| 4     | 32              | 21             | 10      |
