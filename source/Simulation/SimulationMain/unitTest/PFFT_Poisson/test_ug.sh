#!/bin/bash

file=flash.par.test
rm test_ug.out

total[0]=24; iProcs[0]=2; jProcs[0]=4; kProcs[0]=3
total[1]=24; iProcs[1]=4; jProcs[1]=2; kProcs[1]=3
total[2]=24; iProcs[2]=1; jProcs[2]=8; kProcs[2]=3
total[3]=24; iProcs[3]=8; jProcs[3]=1; kProcs[3]=3

for i in 0 1 2 3
do
  valI=${iProcs[$i]}
  valJ=${jProcs[$i]}
  valK=${kProcs[$i]}

  #create a new flash.par file for UG 3D simulations.
  echo "iProcs = $valI
jProcs = $valJ
kProcs = $valK
iGridSize=32
jGridSize=32
kGridSize=24" > $file

  echo "-------------------------------------------" >> test_ug.out
  cat $file >> test_ug.out
  mpirun -np ${total[$i]} ./flash3 -par_file $file 2>&1 >> test_ug.out

done
rm $file
