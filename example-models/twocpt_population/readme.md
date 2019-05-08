- [Build](#org94956bb)
  - [Edit/Add `cmdstan/make/local`](#orgaf7a243)
  - [Build in `cmdstan`](#org9b8c578)
- [Run](#orgce94d85)


<a id="org94956bb"></a>

# Build


<a id="orgaf7a243"></a>

## Edit/Add `cmdstan/make/local`

```sh
TORSTEN_MPI = 1                                         # flag on torsten's MPI solvers
-include $(MATH)make/setup_torsten.mk
CXXFLAGS += $(CXXFLAGS_MPI) -isystem /usr/local/include # path to MPI library's headers
LDFLAGS += $(LDFLAGS_MPI)
CC=mpicxx                                               # name of mpi compilers
CXX=mpicxx
```


<a id="org9b8c578"></a>

## Build in `cmdstan`

```sh
make ../example-models/twocpt_population/twocpt_population
```


<a id="orgce94d85"></a>

# Run

```sh
mpiexec -n 2 twocpt_population sample data file=twocpt_population.data.R init=twocpt_population.init.R
```
