dnl **********************************************************
dnl **********************************************************
dnl CHECKS FOR MPI
dnl **********************************************************
dnl **********************************************************
AC_DEFUN(FAC_CHECK_MPI,
[
AC_MSG_CHECKING(if MPI is wanted)

AC_ARG_WITH(mpi,
 [  --with-mpi=DIR          root directory path of mpi installation [defaults to
                          /usr/local /usr /opt /usr/community]
  --without-mpi          to disable mpi usage completely],
dnl Action if a withvalue is given at command line
 [if test "$withval" != no ; then
   AC_MSG_RESULT(yes)
   MPI_HOME="$withval"
 else
   AC_MSG_RESULT(no)
   MPI_HOME=""
 fi],
dnl Action there was no with(out)-mpi flag on command line
[
AC_MSG_RESULT(yes)
dnl MPI_HOME is set at the beginning based on prototype

if test ! -f "${MPI_HOME}/include/mpi.h"
 then
        for i in /usr/local /usr /opt /usr/community
        do
                if test -f "${i}/include/mpi.h"
                then
                        MPI_HOME=${i}
                else
                        if test -f "${i}/mpi/include/mpi.h"
                        then
                                MPI_HOME=${i}
                        fi
                fi
        done
 fi
])

dnl Locate MPI, if wanted
if test -n "${MPI_HOME}"
then
        MPI_OLD_LDFLAGS=$LDFLAGS
        MPI_OLD_CPPFLAGS=$CPPFLAGS
        MPI_OLD_LIBS=$LIBS
        LDFLAGS="$LDFLAGS -L${MPI_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${MPI_HOME}/include"
        LIBS="${LIBS} -lmpi"
dnl LIBS="${LIBS} -lfmpich -lmpich"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(mpi, MPI_Init, [mpi_cv_libmpi=yes], [mpi_cv_libmpi=no])
dnl AC_HAVE_LIBRARY(fmpich, [LIBS="$LIBS -lfmpich"])
dnl AC_HAVE_LIBRARY(mpich, [LIBS="$LIBS -lmpich"])
        AC_CHECK_HEADER(mpi.h, [mpi_cv_mpi_h=yes], [mpi_cv_mpi_h=no])
        AC_LANG_RESTORE
        if test "$mpi_cv_libmpi" = "yes" -a "$mpi_cv_mpi_h" = "yes"
        then
                dnl If both library and header were found, use them
                AC_MSG_CHECKING(mpi in ${MPI_HOME})
                AC_MSG_RESULT(ok)
                dnl The compiler flags are already set from the check before if
        else
                dnl If either header or library was not found, revert and bomb
                LDFLAGS="${MPI_OLD_LDFLAGS}"
                CPPFLAGS="${MPI_OLD_CPPFLAGS}"
                LIBS="${MPI_OLD_LIBS}"
                AC_MSG_CHECKING(MPI not in ${MPI_HOME})
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid MPI installation with --with-mpi=DIR or disable MPI usage with --without-mpi)
        fi
fi
])

dnl **********************************************************
dnl **********************************************************
dnl CHECKS FOR MPICH - primarily the mpich compilers (mpif90)
dnl **********************************************************
dnl **********************************************************
dnl Does not exist yet. I am thinking about it.
dnl **********************************************************
dnl **********************************************************
dnl CHECKS FOR HDF4
dnl **********************************************************
dnl **********************************************************
AC_DEFUN(FAC_CHECK_HDF4,
[
AC_MSG_CHECKING(which HDF4 installation)
if test ! -f "${HDF4_HOME}/include/hdf.h" ; then
        for i in /usr/local /usr /opt /usr/community
        do
                if test -f "${i}/include/hdf.h"
                then
                        HDF4_HOME=${i}
                else
                        if test -f "${i}/hdf/include/hdf.h"
                        then
                                HDF4_HOME=${i}
                        fi
                fi
        done
fi
dnl Locate HDF, if wanted
if test -n "${HDF4_HOME}"
then
        HDF4_OLD_LDFLAGS=$LDFLAGS
        HDF4_OLD_CPPFLAGS=$CPPFLAGS
        HDF4_OLD_LIBS=$LIBS
        LDFLAGS="$LDFLAGS $HDF4_LDFLAGS"
        CPPFLAGS="$CPPFLAGS $HDF4_CPPFLAGS"
        LIBS="${LIBS} ${HDF4_LIBS}"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(z, inflateEnd, [zlib_cv_libz=yes], [zlib_cv_libz=no])
dnl AC_HAVE_LIBRARY(z, [LIBS="$LIBS -lz"])
dnl AC_HAVE_LIBRARY(jpeg, [LIBS="$LIBS -ljpeg"])
        AC_CHECK_HEADER(hdf.h, [hdf_cv_hdf_h=yes], [zlib_cv_hdf_h=no])
        AC_LANG_RESTORE
        if test "$zlib_cv_libz" = "yes" -a "$hdf_cv_hdf_h" = "yes"
        then
                dnl If both library and header were found, use them
                AC_MSG_CHECKING(hdf4 in ${HDF4_HOME})
                AC_MSG_RESULT(ok)
                dnlThe compiler flags are already set from the check before if
        else
                dnl If either header or library was not found, revert and bomb
                LDFLAGS="${HDF4_OLD_LDFLAGS}"
                CPPFLAGS="${HDF4_OLD_CPPFLAGS}"
                LIBS="${HDF4_OLD_LIBS}"
                AC_MSG_CHECKING(HDF in ${HDF4_HOME} or default test locations)
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid HDF4 installation with --with-hdf=DIR or disable hdf usage with --without-hdf)
        fi
fi
])

dnl **********************************************************
dnl **********************************************************
dnl CHECKS FOR HDF5
dnl **********************************************************
dnl **********************************************************
AC_DEFUN(FAC_CHECK_HDF5,
[
AC_MSG_CHECKING(if hdf5 parallel enabled)
AC_ARG_ENABLE(hdf5parallel,
  [  --enable-hdf5parallel      Uses hdf5 in parallel,
  --disable-hdf5parallel      Does not run in parallel, serial is default],
dnl Enable given at command line. Disable assumes serial
  [
   if test "$enableval" = yes; then
	case $PROTO in
	AIX)
		HDF5_HOME="/usr/local/hdf5/hdf5-1.4.1/parallel"

		CFLAGS="$CFLAGS \$(HDFPATH)/include"
		CFLAGS_DEBUG="$CFLAGS_DEBUG \$(HDFPATH)/include"

		BRENAME="-brename:.h5_open_file_for_read_,.h5_open_file_for_read \
     		    -brename:.h5_read_header_info_,.h5_read_header_info \
       		    -brename:.h5_read_lrefine_,.h5_read_lrefine \
       		    -brename:.h5_read_nodetype_,.h5_read_nodetype \
	            -brename:.h5_read_gid_,.h5_read_gid \
	            -brename:.h5_read_coord_,.h5_read_coord \
	            -brename:.h5_read_size_,.h5_read_size \
	            -brename:.h5_read_bnd_box_,.h5_read_bnd_box \
	            -brename:.h5_read_unknowns_,.h5_read_unknowns \
  	            -brename:.h5_close_file_,.h5_close_file \
	            -brename:.h5_initialize_file_,.h5_initialize_file \
	            -brename:.h5_write_header_info_,.h5_write_header_info \
                    -brename:.h5_write_lrefine_,.h5_write_lrefine \
                    -brename:.h5_write_nodetype_,.h5_write_nodetype \
                    -brename:.h5_write_gid_,.h5_write_gid \
                    -brename:.h5_write_coord_,.h5_write_coord \
                    -brename:.h5_write_size_,.h5_write_size \
                    -brename:.h5_write_bnd_box_,.h5_write_bnd_box \
                    -brename:.h5_write_unknowns_,.h5_write_unknowns \
                    -brename:.h5_write_header_info_sp_,.h5_write_header_info_sp \
                    -brename:.h5_write_coord_sp_,.h5_write_coord_sp \
                    -brename:.h5_write_size_sp_,.h5_write_size_sp \
                    -brename:.h5_write_bnd_box_sp_,.h5_write_bnd_box_sp \
                    -brename:.h5_write_unknowns_sp_,.h5_write_unknowns_sp"

		HDF5_LIBS="-L \$(HDFPATH)/lib -lhdf5 -L /usr/local/lib -lz -bmaxdata:0x80000000"
		AC_MSG_RESULT(yes - config setup for AIX)
	;;
	IRIX64)
		AC_MSG_RESULT(No default hdf5 parallel config for IRIX)
	;;
	Linux)
		AC_MSG_RESULT(No default hdf5 parallel config for Linux)
	;;
	TFLOPS)
		CC=mpicc
		CXX=mpiCC
		F77=mpif90
		FC=mpif90
		LINK=mpicc

		HDF5_HOME="/usr/community/hdf5/hdf5-1_2_1"
		ROMIO_PATH="/usr/community/mpi-io/romio/romio_1.0"

		CFLAGS="$CFLAGS -I \$(HDFPATH) -I $ROMIO_PATH/include"
		CFLAGS_DEBUG="$CFLAGS_DEBUG -I \$(HDFPATH) -I $ROMIO_PATH/include"

		HDF5_LIBS="-L $(HDFPATH)/lib -lhdf5"
		AC_MSG_RESULT(yes - config setup for TFLOPS)
	;;
	UNICOS)
		AC_MSG_RESULT(No hdf5 parallel config UNICOS)
	;;
	*)
	AC_MSG_RESULT(Confused - no prototype for hdf5parallel)
	;;
	esac
  fi
  if test "$enableval" != yes; then
	AC_MSG_RESULT(no)
  fi
],
[
  AC_MSG_RESULT(no)
])

AC_MSG_CHECKING(which HDF5 installation)
dnl HDF5_HOME is set at the beginning based on prototype
if test ! -f "${HDF5_HOME}/include/hdf.h"
 then
        for i in /usr/local /usr /opt /usr/community
        do
                if test -f "${i}/include/hdf.h"
                then
                        HDF5_HOME=${i}
                else
                        if test -f "${i}/hdf/include/hdf.h"
                        then
                                HDF5_HOME=${i}
                        fi
                fi
        done
fi
AC_MSG_RESULT($HDF5_HOME)

dnl Locate HDF5, if wanted
if test -n "${HDF5_HOME}"
then
        HDF5_OLD_LDFLAGS=$LDFLAGS
        HDF5_OLD_CPPFLAGS=$CPPFLAGS
        HDF5_OLD_LIBS=$LIBS
        LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
        CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
        LIBS="${LIBS} ${HDF5_LIBS}"

        AC_LANG_SAVE
        AC_LANG_C

        AC_CHECK_LIB(hdf5, h5_initialize_file, [hdf5_cv_libhdf5=yes], [hdf5_cv_libhdf5=no])
        AC_CHECK_HEADER(hdf5.h, [hdf5_cv_hdf5_h=yes], [hdf5_cv_hdf5_h=no])
        AC_LANG_RESTORE
        if test "$hdf5_cv_libhdf5" = "yes" -a "$hdf5_cv_hdf5_h" = "yes"
        then
                dnl If both library and header were found, use them
                AC_MSG_CHECKING(hdf5 in ${HDF5_HOME})
		HDF_HOME="${HDF5_HOME}"
                AC_MSG_RESULT(ok)
                dnl The compiler flags are already set from the check before if
        else
                dnl If either header or library was not found, revert and bomb
                LDFLAGS="${HDF5_OLD_LDFLAGS}"
                CPPFLAGS="${HDF5_OLD_CPPFLAGS}"
                LIBS="${HDF5_OLD_LIBS}"
                AC_MSG_CHECKING(HDF in ${HDF5_HOME} or default test locations)
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid HDF5 installation with --with-hdf=DIR or disable hdf usage with --without-hdf)
        fi
fi
])

dnl **********************************************************
dnl **********************************************************
dnl Which HDF am I going to use: HDF4 or HDF5
dnl **********************************************************
dnl **********************************************************
AC_DEFUN(FAC_WHICH_HDF,
[
AC_MSG_CHECKING(checking which HDF to use)

dnl Evaluate hdf4 flag
AC_ARG_WITH(hdf4,
[  --with-hdf4=DIR         root directory path of hdf installation [defaults to
                          /usr/local /usr /opt /usr/community]
  --without-hdf4         to disable hdf usage completely],
[if test "$withval" != no ; then
	TEST_HDF4="yes"
	if test "$withval" != yes ; then
		HDF4_HOME="$withval"
	fi
	AC_MSG_RESULT(HDF4 will try  $HDF4_HOME)
	FAC_CHECK_HDF4
else
	TEST_HDF4="no"
fi],
[ TEST_HDF4="no" ])

dnl Evaluate hdf5 flag
AC_ARG_WITH(hdf5,
[  --with-hdf5=DIR         root directory path of hdf installation [defaults to
                          /usr/local /usr /opt /usr/community]
  --without-hdf5          to disable hdf usage completely],
dnl Argument given
[if test "$withval" != no ; then
	if test "$TEST_HDF4" != no ; then
		AC_MSG_ERROR(Trying to use both HDF4 and HDF5)
	else
		if test "$withval" != yes; then
			HDF5_HOME="$withval"
		fi
		TEST_HDF5="yes"
		AC_MSG_RESULT(HDF5 will try $HDF5_HOME)
		FAC_CHECK_HDF5
	fi
else
	TEST_HDF5="no"
fi],
dnl Argument not given
[ TEST_HDF5="no" ])

if test "$TEST_HDF4" = no ; then
	if test "$TEST_HDF5" = no ; then
                AC_MSG_RESULT(no HDF)
        fi
fi

])

