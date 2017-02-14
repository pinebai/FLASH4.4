#!/bin/bash -x
set -e
_INSTALL="${HOME}/software/libNBC/1.1.1"
_PACKAGE="libNBC-1.1.1"
_COMPFLAGS="-O2 -g"
_PACKAGE_URL="http://htor.inf.ethz.ch/research/nbcoll/libnbc/${_PACKAGE}.tar.gz"

patches[0]="libNBC-undef_schedule.patch"
patches[1]="libNBC-function_args.patch"
patches[2]="libNBC-ANSI_C.patch"
patches[3]="libNBC-printf.patch"
patches[4]="libNBC-dcmf_missing_datatype.patch"
patches[5]="libNBC-disable_dcmf.patch"

wget ${_PACKAGE_URL}
tar -zxvf ${_PACKAGE}.tar.gz
cd ${_PACKAGE}

for patch in "${patches[@]}"; do
  patch -p1 -b -i ../"${patch}"
done

MPICXX="mpicc" ; export MPICXX
./configure \
--prefix="${_INSTALL}" \
2>&1 | tee ../${_PACKAGE}_configure.out

make 2>&1 | tee ../${_PACKAGE}_make.out

make install 2>&1 | tee ../${_PACKAGE}_make_install.out

#Keep a record of how I build libNBC.
mkdir ${_INSTALL}/build
for patch in "${patches[@]}"; do
  cp ../"${patch}" ${_INSTALL}/build
done
cp ../"$0" ${_INSTALL}/build
