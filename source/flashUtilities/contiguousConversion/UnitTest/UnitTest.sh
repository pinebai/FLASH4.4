#!/bin/bash

ln -s ../ut_conversionInterface.F90
ln -s ../ut_convertToArrayIndicies.F90
ln -s ../ut_convertToMemoryOffset.F90

make -f UnitTest_Makefile
echo ""
./UnitTest_ContiguousConversion
echo ""
make clean -f UnitTest_Makefile

rm ut_conversionInterface.F90
rm ut_convertToArrayIndicies.F90
rm ut_convertToMemoryOffset.F90
