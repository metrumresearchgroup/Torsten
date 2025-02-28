.PHONY: clean-torsten
clean-torsten:
	@echo '  removing torsten test executables'
	@$(RM) $(call findfiles,stan/math/torsten/test,*_test$(EXE))
	@$(RM) $(call findfiles,stan/math/torsten/test,*_test.d)
	@$(RM) $(call findfiles,stan/math/torsten/test,*_test.d~)
	@$(RM) $(call findfiles,stan/math/torsten/test,*_test.xml)
	@$(RM) $(call findfiles,stan/math/torsten/test,*.o)
