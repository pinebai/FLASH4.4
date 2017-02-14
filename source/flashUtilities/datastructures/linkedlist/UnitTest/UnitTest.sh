#!/bin/bash
export OMP_NUM_THREADS=4
ln -s ../ut_listMethods.includeF90

make -f UnitTest_Makefile
echo ""
./UnitTest_List
#valgrind ./UnitTest_List --leak-check=full -v
echo ""
make clean -f UnitTest_Makefile

rm ut_listMethods.includeF90
