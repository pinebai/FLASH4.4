# Makefile for flashUtilities

flashUtilities += ut_sysMemData.o \
	ut_sysMemInterface.o \
	ut_sysMemStats.o \
	ut_sysMemSummaryStats.o

# The following GNU make special prevents that apparent dependencies
# on the file iso_c_binding.mod, which does usually not actually exist
# in the object directory but refers to the ISO_C_BINDING module known
# to FORTRAN compilers internally, trigger unnecessary recompilation
# of files that refer to the ISO_C_BINDING module.
.SECONDARY: iso_c_binding.mod

