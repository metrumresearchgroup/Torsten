# By default we use DQ, but user can choose AD for Jacobian
# computation in CVODES' newton iterations.
ifdef TORSTEN_CVS_JAC_AD
   CXXFLAGS += -DTORSTEN_CVS_JAC_AD
endif

# turning on braid implies using MPI
ifdef TORSTEN_BRAID
  TORSTEN_MPI = 2
  CXXFLAGS_MPI += -DTORSTEN_BRAID -isystem $(BRAID_PATH)/braid
  LDLIBS_MPI += $(LIBBRAID)
  LDLIBS_MPI += $(BRAID_PATH)/braid/libbraid.a
endif

# Torsten's MPI may have different setup than Stan's MPI
ifdef TORSTEN_MPI
  # LIBMPI ?=
  CXXFLAGS_MPI += -DTORSTEN_MPI
  ifeq ($(TORSTEN_MPI), 2)
    CXXFLAGS_MPI += -DTORSTEN_MPI_DYN # dynamic load balance
  endif
  CXXFLAGS += $(CXXFLAGS_MPI)
  LDFLAGS += $(LDFLAGS_MPI)
  CC=mpicxx
  CXX=mpicxx
  # LDFLAGS_MPI ?=
endif

# By default Torsten solvers use Torsten's ODE integrator,
# but user can choose using Stan's integrators.
ifdef TORSTEN_USE_STAN_ODE
  CXXFLAGS += -DTORSTEN_USE_STAN_ODE
endif

