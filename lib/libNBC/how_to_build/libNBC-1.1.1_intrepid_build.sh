#!/bin/bash -x
set -e
_TRANSPORT="mpi" # "dcmf" # "mpi"
_INSTALL="/intrepid-fs0/users/cdaley/persistent/software/libNBC/1.1.1/xl/${_TRANSPORT}"
_PACKAGE="libNBC-1.1.1"
_COMPFLAGS="-O3 -qnohot -g"

patches[0]="libNBC-undef_schedule.patch"
patches[1]="libNBC-function_args.patch"
patches[2]="libNBC-ANSI_C.patch"
patches[3]="libNBC-printf.patch"
patches[4]="libNBC-dcmf_missing_datatype.patch"
if [ "${_TRANSPORT}" != "dcmf" ]; then
    patches[5]="libNBC-disable_dcmf.patch"
fi


tar -zxvf ${_PACKAGE}.tar.gz
cd ${_PACKAGE}

for patch in "${patches[@]}"; do
  patch -p1 -b -i ../"${patch}"
done

MPICXX="mpixlc_r" ; export MPICXX
./configure --host=powerpc-bgp-linux \
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
