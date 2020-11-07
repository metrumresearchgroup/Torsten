.PHONY: doxygen
doxygen:
	mkdir -p doc/api
	doxygen doxygen/doxygen.cfg

.PHONY: clean clean-doxygen clean-deps clean-all
# clean:
# 	@echo '  removing torsten test executables'
# 	@$(RM) $(call findfiles,stan/math/torsten/test,*_test$(EXE))
# 	@$(RM) $(call findfiles,stan/math/torsten/test,*_test.d)
# 	@$(RM) $(call findfiles,stan/math/torsten/test,*_test.d.*)
# 	@$(RM) $(call findfiles,stan/math/torsten/test,*_test.xml)
# 	@$(RM) $(call findfiles,stan/math/torsten/test,*.o)

clean-doxygen:
	@echo '  removing doxygen'
	$(RM) -r doc/api

# clean-deps:
# 	@echo '  removing dependency files'
# 	@$(RM) $(call findfiles,stan,*.d)
# 	@$(RM) $(call findfiles,test,*.d)
# 	@$(RM) $(call findfiles,lib,*.d)
# 	@$(RM) $(call findfiles,stan,*.d.*)
# 	@$(RM) $(call findfiles,test,*.d.*)
# 	@$(RM) $(call findfiles,lib,*.d.*)
# 	@$(RM) $(call findfiles,stan,*.dSYM)

clean-all: clean clean-doxygen clean-deps clean-libraries
