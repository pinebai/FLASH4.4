# GNU make makefile.
# Builds unit test with one of several compilers:
# gnu, intel, ibm, absoft, lahey and sun.  The default is gnu.
# If you would like to select a specific compiler choose:
# 1). make gnu
# 2). make intel
# 3). make ibm
# 4). make absoft
# 5). make lahey
# 6). make sun

SRC     = Driver.F90 \
	PersonNode.F90 \
	PersonList.F90 \
	GroupNode.F90 \
	GroupList.F90 \
	UnitTest.F90
OBJ     = $(SRC:.F90=.o)
EXE      = UnitTest_ListOfList

#DEFAULT: GNU compiler.
gnu:	FC = gfortran
gnu:	FFLAGS = -fdefault-real-8 -fdefault-double-8 -pedantic \
	-Wall -Wextra -Wconversion -Wunderflow \
	-ffpe-trap=invalid,zero,overflow -fbounds-check -g
gnu:	all


#Intel compiler.
intel:	FC = ifort
intel:	FFLAGS = -r8 -i4 -check bounds -check format \
	-check output_conversion -warn all -real_size 64 \
	-check uninit -traceback -fp-stack-check -g
intel:	all


#IBM compiler.
ibm:	FC = xlf2003
ibm:	FFLAGS = -qfloat=rngchk -qcheck -qrealsize=8 -qnosave \
	-qflttrap=enable:invalid:nanq:overflow:zerodivide -g
ibm:	all


#Absoft compiler.
absoft:	FC = absf90
absoft:	FFLAGS = -g -N113 -en -m64 -et \
	-trap=INVALID,DIVBYZERO,OVERFLOW -B111 -Rb -Rc -Rs -safefp +B71
absoft:	all


#Lahey compiler
lahey:	FC = lf95
lahey:	FFLAGS = -g --trace --trap --chk[aes] -CcdRR8
lahey:	all


#Sun Studio compiler.
sun:	FC = sunf95
sun:	FFLAGS = -xcheck=%all -ftrap=invalid,overflow,division \
	-xtypemap=real:64,double:64 -g
sun:	all


all:    $(OBJ)
	$(FC) $(FFLAGS) -o $(EXE) $(OBJ)

%.o:    %.F90
	$(FC) -c $< $(FFLAGS) -o $@


PersonList.o: ut_listMethods.includeF90
GroupList.o: ut_listMethods.includeF90


clean:
	rm -f $(OBJ) $(EXE) *.mod
