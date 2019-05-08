- [Build](#orgba67c47)
  - [Edit/Add `cmdstan/make/local`](#org8585aae)
  - [Build in `cmdstan`](#org532eb39)
- [Run](#org7bf6c43)


<a id="orgba67c47"></a>

# Build


<a id="org8585aae"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org532eb39"></a>

## Build in `cmdstan`

```sh
make ../example-models/lotka_volterra_ode_group_model/lv_group
```


<a id="org7bf6c43"></a>

# Run

```sh
mpiexec -n 2 lv_group sample data file=lv_group.data.R
```
