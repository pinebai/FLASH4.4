#!/bin/bash

SEARCH_DIRS="/lib /usr"
SEARCH_LIBS="mpi hdf5 fftw netcdf papi"

function findlib()
{
   LIBNAME=$1
   for DIR in $SEARCH_DIRS; do
       for m in `find $DIR -type d -maxdepth 5 -name lib -ipath *${LIBNAME}* 2>/dev/null`; do 
           echo -n "   ";
           dirname $m;
       done;
   done;
}

echo "Searching for the following libraries: $SEARCH_LIBS"
for lib in $SEARCH_LIBS; do
    echo "$lib found at:"
    findlib $lib;
done
