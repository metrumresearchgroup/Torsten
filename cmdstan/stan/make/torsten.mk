# Adding Torsten functions making MPL list too long, need adjust list size
CXXFLAGS += -DBOOST_MPL_CFG_NO_PREPROCESSED_HEADERS -DBOOST_MPL_LIMIT_LIST_SIZE=30

ifeq ($(STANC2),)
TORSTEN_STANC3_VERSION=torsten_v0.91.1

# supply path to STANC3 to run custom stanc3 transpiler
#  STANC3=../stanc3
endif
