#!/bin/bash

ln -s ../iso.c .
ln -s ../cell_table.h .
ln -s ../ut_contourSurfaceAreaBlock.F90 .
ln -s ../ut_contourSurfaceInterface.F90 .
cp ../../general/mangle_names.h .


if make -f UnitTest_Makefile
then
  # leave faild builds behind
  echo
  ./UnitTest_contourSurface
  retval=$?

  make clean -f UnitTest_Makefile

  rm iso.c
  rm cell_table.h
  rm ut_contourSurfaceAreaBlock.F90
  rm ut_contourSurfaceInterface.F90
  rm mangle_names.h

fi
exit $retval

