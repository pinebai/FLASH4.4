#!/bin/bash
ln -s ../ut_listMethods.includeF90

make
echo ""
./UnitTest_ListOfList
echo ""
make clean

rm ut_listMethods.includeF90
