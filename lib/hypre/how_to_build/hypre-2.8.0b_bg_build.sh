#!/bin/bash
set -e

# This script builds four different versions of HYPRE on Intrepid
# (BG/P) and Mira (BG/Q).  It will create MPI and MPI+OpenMP versions
# of HYPRE built against 1). internal LAPACK and BLAS libraries and
# 2). IBMs tuned ESSL library:
# http://www-03.ibm.com/systems/software/essl/

# Resources:
# https://computing.llnl.gov/LCdocs/math/math.pdf
# https://www.alcf.anl.gov/user-guides/bgp-faqs-compiling-and-linking
# https://code.google.com/p/blopex/wiki/HypreInstallLinux
# http://bluegene.bnl.gov/comp/bgp/libraries.shtml
# https://wiki.fysik.dtu.dk/gpaw/install/BGP/building_with_gcc_on_surveyor.html
# https://github.com/ahmadia/elemental/blob/master/cmake/toolchains/BGP-Shaheen-xl.cmake
# https://www.ccs.uky.edu/Events/IBM_Workshop_07/workshop_files/HPC_session_4.pdf
# https://code.google.com/p/m-a-d-n-e-s-s/wiki/UsingALCF

# Notes:
#
# 1. Patch the configure script so that the Blue Gene ESSL library,
# which is named libesslbg.a and not libessl.a, is recognized as an
# ESSL library.  This ensures the macro HYPRE_USING_ESSL is defined.
#
# 2. Manually add CPPFLAGS to CFLAGS and CXXFLAGS because CPPFLAGS is
# ignored during the compilation of HYPRE.  This ensures the ESSL
# header file is found.
#
# 3. Add the libraries "xlf90_r xlfmath xlomp_ser" to the
# --with-blas-libs configure option so that we can build the HYPRE
# test programs.  The HYPRE test programs are written in C and are
# built by running "make test".  These additional libraries are not
# needed when building FLASH against HYPRE because FLASH is linked
# with a Fortran compiler.

_PACKAGE="hypre-2.8.0b"
_MAKE_JOBS=2
# config.guess returns the same build environment on Mira and Intrepid login nodes
# Intrepid login node processors: PowerPC 970MP (64-bit)
# Mira login node processors: Power 740 (64-bit)
_BUILD="powerpc64-unknown-linux-gnu"

name="$(hostname)"
if [[ "$name" == intrepid* ]] || [[ "$name" == challenger* ]]; then
  echo "BG/P platform!"
  _HOST="powerpc-bgp-linux"
  _ARCH_OPTION="-qarch=450"
  _INSTALL="/intrepid-fs0/users/cdaley/persistent/software/hypre/2.8.0b/xl"
  _ESSL_PATH="/soft/apps/current/ESSL-4.4.1-1"
  _ESSL_LIBS="${_ESSL_PATH}/lib ${IBM_MAIN_DIR}/xlf/bg/11.1/lib"
elif [[ "$name" == mira* ]] || [[ "$name" == cetus* ]]; then
  echo "BG/Q platform!"
  _HOST="powerpc64-bgq-linux"
  _ARCH_OPTION="-qarch=qp"
  _INSTALL="/gpfs/mira-fs0/projects/TurbNuclComb_esp/cdaley/software/V1R2M0_efix38/hypre/2.8.0b/xl"
  _ESSL_PATH="/soft/libraries/essl/5.1.1-0"
  _ESSL_LIBS="${_ESSL_PATH}/lib64 ${IBM_MAIN_DIR}/xlf/bg/14.1/lib64"
else
  echo "Platform not recognized. Exiting"
  exit -1
fi


# Shouldn't need to edit below this point

tar -zxvf ${_PACKAGE}.tar.gz
cd ${_PACKAGE}/src

patches[0]="hypre-2.8.0b-esslbg.patch"
for patch in "${patches[@]}"; do
  patch -p2 -b -i ../../"${patch}"
done


export CXX=mpixlcxx_r
export CC=mpixlc_r
export F77=mpixlf77_r


# *** Build an MPI-only version of HYPRE ***
_COMPFLAGS="-g -O2 ${_ARCH_OPTION}"
export FFLAGS="${_COMPFLAGS}"
export F77FLAGS="${_COMPFLAGS}"

# Build against ESSL
export CPPFLAGS="-I${_ESSL_PATH}/include"
export CXXFLAGS="${CPPFLAGS} ${_COMPFLAGS}"
export CFLAGS="${CPPFLAGS} ${_COMPFLAGS}"
_PREFIX="${_INSTALL}/essl"
./configure \
    --build="${_BUILD}" --host="${_HOST}" \
    --with-blas-libs="esslbg xlf90_r xlfmath xlomp_ser" \
    --with-blas-lib-dirs="${_ESSL_LIBS}" \
    --with-lapack-libs="" \
    --with-lapack-lib-dirs="" \
    --prefix="${_PREFIX}" \
    --with-MPI \
    --with-timing \
    --with-print-errors
make -j ${_MAKE_JOBS} && make -j ${_MAKE_JOBS} test && make install
cp config.log ${_PREFIX}
make distclean

# Build against the internal LAPACK and BLAS libraries
unset CPPFLAGS
export CXXFLAGS="${_COMPFLAGS}"
export CFLAGS="${_COMPFLAGS}"
_PREFIX="${_INSTALL}/default"
./configure \
    --build="${_BUILD}" --host="${_HOST}" \
    --prefix="${_PREFIX}" \
    --with-MPI \
    --with-timing \
    --with-print-errors
make -j ${_MAKE_JOBS} && make -j ${_MAKE_JOBS} test && make install
cp config.log ${_PREFIX}
make distclean


# *** Build an MPI+OpenMP version of HYPRE ***
_OPENMP="-qsmp=omp:noauto"
_COMPFLAGS="${_OPENMP} -g -O2 ${_ARCH_OPTION}"
export FFLAGS="${_COMPFLAGS}"
export F77FLAGS="${_COMPFLAGS}"
export LDFLAGS="${_OPENMP}"

# Build against ESSL
export CPPFLAGS="-I${_ESSL_PATH}/include"
export CXXFLAGS="${CPPFLAGS} ${_COMPFLAGS}"
export CFLAGS="${CPPFLAGS} ${_COMPFLAGS}"
_PREFIX="${_INSTALL}/essl-omp"
./configure \
    --build="${_BUILD}" --host="${_HOST}" \
    --with-openmp \
    --with-blas-libs="esslbg xlf90_r xlfmath" \
    --with-blas-lib-dirs="${_ESSL_LIBS}" \
    --with-lapack-libs="" \
    --with-lapack-lib-dirs="" \
    --prefix="${_PREFIX}" \
    --with-MPI \
    --with-timing \
    --with-print-errors
make -j ${_MAKE_JOBS} && make -j ${_MAKE_JOBS} test && make install
cp config.log ${_PREFIX}
make distclean

# Build against the internal LAPACK and BLAS libraries
unset CPPFLAGS
export CXXFLAGS="${_COMPFLAGS}"
export CFLAGS="${_COMPFLAGS}"
_PREFIX="${_INSTALL}/default-omp"
./configure \
    --build="${_BUILD}" --host="${_HOST}" \
    --with-openmp \
    --prefix="${_PREFIX}" \
    --with-MPI \
    --with-timing \
    --with-print-errors
make -j ${_MAKE_JOBS} && make -j ${_MAKE_JOBS} test && make install
cp config.log ${_PREFIX}
make distclean
