- [Build](#orge9524af)
  - [Edit/Add `cmdstan/make/local`](#org8306633)
  - [Build in `cmdstan`](#org7bd3108)
- [Run](#orgf025eab)
- [Results](#orge4e2a4f)


<a id="orge9524af"></a>

# Build


<a id="org8306633"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org7bd3108"></a>

## Build in `cmdstan`

```sh
make ../example-models/lotka_volterra_ode_group_model/lv_group
```


<a id="orgf025eab"></a>

# Run

```sh
mpiexec -n 2 lv_group sample data file=lv_group.data.R
```


<a id="orge4e2a4f"></a>

# Results

Three chains are run using

-   sequential run using Stan's `rk45` integrator(output `stan_sample.1-4.csv`)
-   sequential run using Torsten's `rk45` integrator(output `pmx_sample.1-4.csv`)
-   MPI run using Torsten's `rk45` group integrator(output `mpi_sample.1-4.csv`)

The wall time of sequential runs and MPI runs(in seconds):

| chain | sequential stan | sequential pmx | MPI pmx |
|----- |--------------- |-------------- |------- |
