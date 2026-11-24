N_TESTS ?= 100
# override this in make/local. If <= 0, N_TESTS + 1 is interpreted as the number of batches
# to group the probability tests into

##
# Any targets in test/ (.d, .o, executable) needs the GTEST flags
##

stan/math/torsten/test/% : CXXFLAGS += $(CXXFLAGS_GTEST)
stan/math/torsten/test/% : CPPFLAGS += $(CPPFLAGS_GTEST)
stan/math/torsten/test/% : INC += $(INC_GTEST)

stan/math/torsten/test/%$(EXE) : stan/math/torsten/test/%.o $(GTEST)/src/gtest_main.cc $(GTEST)/src/gtest-all.o $(MPI_TARGETS) $(TBB_TARGETS)
	$(LINK.cpp) $^ $(LDLIBS) $(OUTPUT_OPTION)

##
# Include dependency files for tests
##
ifneq ($(filter stan/math/torsten/test/%,$(MAKECMDGOALS)),)
-include $(patsubst %$(EXE),%.d,$(filter stan/math/torsten/test/%,$(MAKECMDGOALS)))
endif

############################################################
#
# torsten tests
##

TORSTEN_TESTS := $(subst .cpp,$(EXE),$(shell find stan/math/torsten/test/unit -name *_test.cpp))
$(TORSTEN_TESTS) : $(SUNDIALS_TARGETS)

TORSTEN_FULL_TESTS := $(subst .cpp,$(EXE),$(shell find stan/math/torsten/test/full_unit -name *_test.cpp))
$(TORSTEN_FULL_TESTS) : $(SUNDIALS_TARGETS)

TORSTEN_PERF_TESTS := $(subst .cpp,$(EXE),$(shell find stan/math/torsten/test/performance -name *_test.cpp))
$(TORSTEN_PERF_TESTS) : $(SUNDIALS_TARGETS)

TORSTEN_MPI_DYN_TESTS_DEP := $(subst .cpp,.d,$(shell find stan/math/torsten/test/unit -name mpi_dynamic_load*_test.cpp))
$(TORSTEN_MPI_DYN_TESTS_DEP) : CXXFLAGS_MPI += -DTORSTEN_MPI_DYN

ifdef TORSTEN_MPI
  TORSTEN_MPI_DYN_TESTS := $(subst .cpp,$(EXE),$(shell find stan/math/torsten/test/unit -name mpi_dynamic_load*_test.cpp))
  $(TORSTEN_MPI_DYN_TESTS) : CXXFLAGS_MPI += -DTORSTEN_MPI_DYN
endif

ifneq ($(call ifdef_any_of,TORSTEN_MPI STAN_LANG_MPI),)
#  MPI_TESTS := $(subst .cpp,$(EXE),$(shell find test -name *mpi_*test.cpp))
  GTEST_CXXFLAGS += $(CXXFLAGS_MPI)
endif

############################################################
#
# Target to verify header files within Stan has
# enough include calls
##
HEADER_TESTS := $(addsuffix -test,$(call findfiles,stan,*.hpp))

ifeq ($(OS),Windows_NT)
  DEV_NULL = nul
else
  DEV_NULL = /dev/null
endif

%.hpp-test : %.hpp stan/math/torsten/test/dummy.cpp
	$(COMPILE.cpp) $(CXXFLAGS) -O0 -include $^ -o $(DEV_NULL)

stan/math/torsten/test/dummy.cpp:
	@mkdir -p test
	@touch $@
	@echo "int main() {return 0;}" >> $@

.PHONY: test-headers
test-headers: $(HEADER_TESTS)
