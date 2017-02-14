# name of the library -- make sure it is in the object directory
LIBNAME = ../object/libma28.a

# what the makefile should do -- this is the default thing to do
.PHONY: default
default: $(LIBNAME)

# Import site specific information.
# At setup time, this file will exist in the local library directory.
include  ../Makefile.h

# BUILDFLAG can be passed to the Makefile to choose specific
# compilation flags (one of OPT,TEST,DEBUG).
ifdef BUILDFLAG
	FFLAGS = $(FFLAGS_$(BUILDFLAG))
else
	FFLAGS = $(FFLAGS_OPT)
endif

# command for generating a library
MKLIB = ar -r

# name of the files in the source directory that must be compiled
OBJS = Ma28.o

# dependencies:  the library requires the object files
$(LIBNAME): $(OBJS) 
	$(MKLIB) $(LIBNAME) $(OBJS)

# Dependencies of object files upon source code
# General rule for compiling object files from fortran files
%.o:  %.F90
	$(FCOMP) $(FFLAGS) $(F90FLAGS) -c $<

.PHONY: clean
clean : 
	rm -f $(LIBNAME) $(OBJS)
