#!/bin/bash
set -e
_INSTALL="${HOME}/software/hdf5/1.8.10/xl"
_PACKAGE="hdf5-1.8.10"

tar -zxvf ${_PACKAGE}.tar.gz
cd ${_PACKAGE}

./configure \
CC=mpixlc_r \
FC=mpixlf90_r \
CFLAGS="-O3 -qnohot -g -DIBMR2Fortran -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_AIX -UH5_HAVE_GETPWUID" \
FCFLAGS="-O3 -qnohot -g -WF,-DIBMR2Fortran -WF,-D_LARGEFILE64_SOURCE -WF,-D_FILE_OFFSET_BITS=64 -WF,-D_AIX -WF,-UH5_HAVE_GETPWUID" \
--prefix="${_INSTALL}" \
--enable-parallel \
--enable-production \
--enable-debug=all \
--enable-using-memchecker \
--enable-fortran \
2>&1 | tee ../${_PACKAGE}_configure.out

make 2>&1 | tee ../${_PACKAGE}_make.out

make install 2>&1 | tee ../${_PACKAGE}_make_install.out
