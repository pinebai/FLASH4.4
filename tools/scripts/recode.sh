#!/bin/bash

pushd ~/FLASH3/source/Grid/GridMain/paramesh/Paramesh3
rm -rf *
svn update
rm -f PM3_package/source/amr_*.F.dan 2> /dev/null


# Find all .F files and transform them one at a time
for m in `find . -name '*.F'`; do 
    if ~/FLASH3/tools/scripts/cleanPM.py --input=$m --output=${m}90 ; then
       rm $m ;
       echo Transforming $m ;
    else
       echo *** Unable to transform $m ;
    fi
done
popd

