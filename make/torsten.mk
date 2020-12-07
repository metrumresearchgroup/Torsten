CXXFLAGS += -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_MPL_LIMIT_LIST_SIZE=30
-include $(MATH)make/torsten_setup.mk

ifdef MPI_ADAPTED_WARMUP

CROSS_CHAIN_BOOST_TARGETS ?= $(addsuffix $(LIBRARY_SUFFIX),$(BOOST)/stage/lib/libboost_filesystem $(BOOST)/stage/lib/libboost_regex)
BOOST_LIBRARY_ABSOLUTE_PATH=$(abspath $(BOOST)/stage/lib)
LDFLAGS_CROSS_CHAIN ?= -Wl,-L,"$(BOOST_LIBRARY_ABSOLUTE_PATH)" -Wl,-rpath,"$(BOOST_LIBRARY_ABSOLUTE_PATH)"

CXXFLAGS += -DSTAN_LANG_MPI -DMPI_ADAPTED_WARMUP
LDFLAGS += $(LDFLAGS_CROSS_CHAIN)
CC=mpicxx
CXX=mpicxx

##
# Boost build options
BOOST_PARALLEL_JOBS ?= 2


$(BOOST)/user-config.jam:
	echo "# using boost filesystem & regex:" >> $(BOOST)/user-config.jam

$(CROSS_CHAIN_BOOST_TARGETS): $(BOOST)/user-config.jam
	@mkdir -p $(dir $@)
	cd $(BOOST); ./bootstrap.sh
	cd $(BOOST); ./b2 --user-config=user-config.jam --layout=system --with-filesystem --with-regex -j$(BOOST_PARALLEL_JOBS) variant=release link=shared threading=multi runtime-link=shared hardcode-dll-paths=true dll-path="$(BOOST_LIBRARY_ABSOLUTE_PATH)" cxxstd=11
ifeq ($(OS),Darwin)
	install_name_tool -id @rpath/libboost_filesystem.dylib "$(BOOST)/stage/lib/libboost_filesystem.dylib"
	install_name_tool -id @rpath/libboost_regex.dylib "$(BOOST)/stage/lib/libboost_regex.dylib"
endif

clean-cross-chain-boost-libs:
	@echo '  cleaning boost targets for cross-chain warmup'
	$(RM) $(wildcard $(CROSS_CHAIN_BOOST_TARGETS))
	$(RM) -r $(wildcard $(BOOST)/stage/lib $(BOOST)/bin.v2 $(BOOST)/tools/build/src/engine/bootstrap/ $(BOOST)/tools/build/src/engine/bin.* $(BOOST)/project-config.jam* $(BOOST)/b2 $(BOOST)/bjam $(BOOST)/bootstrap.log)

else

CROSS_CHAIN_BOOST_TARGETS ?=

clean-cross-chain-boost-libs:

endif

.PHONY: build-cross-chain-deps
build-cross-chain-deps: $(CROSS_CHAIN_BOOST_TARGETS)
	@echo '$(CROSS_CHAIN_BOOST_TARGETS)'
	@echo '--- boost cross-chain warmup bindings built ---'
