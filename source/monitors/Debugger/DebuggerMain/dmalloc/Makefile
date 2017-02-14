# The +dmalloc and +ddt setup shortcuts add DMALLOC and DDT macros to
# CDEFINES.
#
# The Makefile fragment below modifies CFLAGS and LIB if CDEFINES
# contains DMALLOC or DDT.  Manipulation of CFLAGS and LIB in this way
# is extremely ugly, but it ensures that the malloc replacement
# library is last in the link line.  This cannot be guaranteed if
# using the library capabililty of FLASH.

ifneq (,$(findstring DMALLOC,$(CDEFINES)))
    CFLAGS += -I$(DMALLOC_PATH)/include
    LIB += -L$(DMALLOC_PATH)/lib -ldmalloc
endif

ifneq (,$(findstring DDT,$(CDEFINES)))
    LIB += -L$(LIB_DDT) -ldmalloc
endif

# Memory Debugging with DDT
# (http://www.alcf.anl.gov/book/export/html/10470)
#
# In order to perform memory debugging with DDT on Blue Gene systems,
# you must statically link your program with DDT's dmalloc library.
# On Intrepid, Challenger, and Surveyor, the library for the default
# version of ddt is installed at /soft/apps/fen/ddt/lib/32.  To use
# it, include the following in your link line:
#
#    -L/soft/apps/fen/ddt/lib/32 -ldmalloc
#
# If you are using the standard compiler wrappers (e.g. mpixlc or
# mpicc), the "-ldmalloc" should be at the end of your link line.  If
# you are listing system libraries explicitly, make sure that
# -ldmalloc comes before -lc.
#
# Note: If you are not using the default version of ddt, please use
# the appropriate path from /soft/apps/fen/ddt-* instead of
# /soft/apps/fen/ddt
#
# During linking, you may encounter errors due to multiple definitions
# of some symbols such as malloc, free, and realloc.  If you have
# confirmed your link order is correct, these errors may be bypassed
# by using the linker option "-z muldefs" or the XL compiler option
# "-Wl,--allow-multiple-definition".
#
# After successfully linking your program with dmalloc, start your
# debugging session as described above but be sure to check the box
# labelled "Memory Debugging".
