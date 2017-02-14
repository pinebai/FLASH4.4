# GNU make makefile.
# Builds pointer passing program with one of several compilers:
# gnu, intel, ibm and sun.  The default is gnu.  If you would like 
# to select a specific compiler choose:
# 1). make gnu
# 2). make intel
# 3). make ibm
# 4). make sun

CXXSRC   = chombo_interface.cpp language_interface.cpp
CXXOBJ   = $(CXXSRC:.cpp=.o)
FSRC     = flash_subroutines.F90 main.F90
FOBJ     = $(FSRC:.F90=.o)
EXE      = unit_test_1


#DEFAULT: GNU compiler.
#Success with gfortran 4.4.2 and 4.5.0 (14th Jan 2010 snapshot)
#GCC bug report 40962 indicates we need gfortran > 4.4.0.
gnu:	FC = gfortran
gnu:	FFLAGS = -fdefault-real-8 -fdefault-double-8 -pedantic \
	-Wall -Wextra -Wconversion -Wunderflow \
	-ffpe-trap=invalid,zero,overflow -fbounds-check -g
gnu:	CXX = g++
gnu:	CXXFLAGS = -pedantic -Wall -Wextra -g
gnu:	LDFLAGS = -lstdc++
gnu:	all


#Intel compiler.
#Success with Intel 10.1.
intel:	FC = ifort
intel:	FFLAGS = -r8 -i4 -check bounds -check format \
	-check output_conversion -warn all -real_size 64 \
	-check uninit -traceback -fp-stack-check -g
intel:	CXX = icc
intel:	CXXFLAGS = -debug extended -g
intel:	LDFLAGS = -static -lstdc++
intel:	all


#IBM compiler.
#Success with xlf 11.1 and xlc++ 9.0 on intrepid.alcf.anl.gov (BG/P).
#The LDFLAGS flags need to be updated whenever there are compiler
#upgrades.  These flags were correct on 28 Apr 2011.
ibm:	FC = xlf2003_r
ibm:	FFLAGS = -qfloat=rngchk -qcheck -qrealsize=8 -qnosave \
	-qflttrap=enable:invalid:nanq:overflow:zerodivide -g
ibm:	CXX = xlc++_r
ibm:	CXXFLAGS = -qfloat=rngchk -qcheck=all \
	-qflttrap=enable:invalid:nanq:overflow:zerodivide -g
ibm:	LDFLAGS = -Wl,--static \
	-L/soft/apps/ibmcmp-dec2010/vacpp/bg/9.0/lib -libmc++ -lstdc++
ibm:	all


#Sun Studio compiler.
#Success with sun studio 12u1.
sun:	FC = sunf95
sun:	FFLAGS = -xcheck=%all -ftrap=invalid,overflow,division \
	-xtypemap=real:64,double:64 -g
sun:	CXX = sunCC
sun:	CXXFLAGS = -xcheck=%all -ftrap=invalid,overflow,division -g
sun:	LDFLAGS = -lCstd
sun:	all


all:    $(CXXOBJ) $(FOBJ)
	$(FC) $(FFLAGS) -o $(EXE) $(FOBJ) $(CXXOBJ) $(LDFLAGS)

#Static Pattern Rules.
$(CXXOBJ): %.o: %.cpp
	$(CXX) -c $< $(CXXFLAGS) -o $@

$(FOBJ): %.o: %.F90
	$(FC) -c $< $(FFLAGS) -o $@

chombo_interface.o: chombo_interface.hpp
language_interface.o: Flash.h
flash_subroutines.o: Flash.h


clean:
	rm -f $(FOBJ) $(CXXOBJ) $(EXE) *.mod *__genmod.*
