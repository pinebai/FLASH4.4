#!/bin/bash

#Just runs PFFT simulation many times with the same flash.par 
#file but various numbers of processors.  In the simulation we are doing 
#a forward and back FFT, and checking results against original 
#values.  We will hopefully get very small values from the simulation.

#Using various numbers of processors stresses the mapping 
#from FLASH grid to PFFT grid, and back again.

parfile=flash.par.test
outfile=test_pm.out
rm $outfile

#lrefine_min = lrefine_max = 1
echo "lrefine_min = 1
lrefine_max = 1
convertToConsvdInMeshInterp           = .false." > $parfile

for i in 1
do
  echo "---------------------------------------------------" >> $outfile
  cat $parfile >> $outfile
  mpirun -np $i ./flash3 -par_file $parfile 2>&1 >> $outfile
done

echo "---------------------------------------------------" >> $outfile
echo "---------------------------------------------------" >> $outfile
echo "---------------------------------------------------" >> $outfile


#lrefine_min = lrefine_max = 2
echo "lrefine_min = 2
lrefine_max = 2
convertToConsvdInMeshInterp           = .false." > $parfile

for i in 1 2 3 4 5 6 7 8
do
  echo "---------------------------------------------------" >> $outfile
  cat $parfile >> $outfile
  mpirun -np $i ./flash3 -par_file $parfile 2>&1 >> $outfile
done


#lrefine_min = lrefine_max = 3
echo "---------------------------------------------------" >> $outfile
echo "---------------------------------------------------" >> $outfile
echo "---------------------------------------------------" >> $outfile
echo "lrefine_min = 3
lrefine_max = 3
convertToConsvdInMeshInterp           = .false." > $parfile

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
do
  echo "---------------------------------------------------" >> $outfile
  cat $parfile >> $outfile
  mpirun -np $i ./flash3 -par_file $parfile 2>&1 >> $outfile
done


grep "the result is" test_pm.out > tmp.1
echo "Line containing max values in column 6 is:"
sort -g +6 -r tmp.1 | head -1
echo "Line containing max values in column 7 is:"
sort -g +7 -r tmp.1 | head -1
echo " --- Very small values, e.g. less than 1E-13, indicate success ---"


rm tmp.1
rm $parfile
