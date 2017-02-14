Grid += flash_ftypes.o chombo_f_c_api.o chombo_f_c_interface.o chombo_uniform_grid.o chombo_adaptive_grid.o AMRLevelFlash.o AMRLevelFlashFactory.o

#It seems like the FLASH setup script does not figure out
#header file dependencies.
AMRLevelFlash.o : AMRLevelFlash.C AMRLevelFlash.h
AMRLevelFlashFactory.o : AMRLevelFlashFactory.C AMRLevelFlashFactory.h
AMRLevelFlashFactory.h : AMRLevelFlash.h
chombo_adaptive_grid.o : chombo_adaptive_grid.C chombo_adaptive_grid.h flash_ctypes.h
chombo_uniform_grid.o: chombo_uniform_grid.C chombo_uniform_grid.h flash_ctypes.h
chombo_f_c_api.o : chombo_f_c_api.C chombo_f_c_api.h
