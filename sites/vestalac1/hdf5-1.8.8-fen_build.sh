#!/bin/bash
_INSTALL="/veas-fs0/cdaley/software/hdf5-fen/1.8.8/xl"
_PACKAGE="hdf5-1.8.8"
_COMPFLAGS="-O3 -qnohot -g -DIBMR2Fortran -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_AIX -UH5_HAVE_GETPWUID"

tar -zxvf ${_PACKAGE}.tar.gz
cd ${_PACKAGE}

./configure \
CC=/soft/compilers/ibmcmp-feb2012/vac/bg/12.1/bin/xlc_r \
CFLAGS="${_COMPFLAGS}" \
--prefix="${_INSTALL}" \
--enable-production \
--enable-debug=all \
--enable-using-memchecker \
2>&1 | tee ../${_PACKAGE}_configure.out

make 2>&1 | tee ../${_PACKAGE}_make.out

make install 2>&1 | tee ../${_PACKAGE}_make_install.out
