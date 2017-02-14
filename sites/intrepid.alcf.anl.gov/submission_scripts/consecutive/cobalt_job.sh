#!/bin/bash -x
# Asssumptions
#  * the qsub line requests a script mode, e.g.
#    qsub ... --mode script cobalt_job.sh
#  * the data files and flash4 binary are in the current working directory
#  * the thread experiment directories do not exist

FLASH_BINARY='flash4'
NODES='256'
THREAD_EXPERIMENTS='4 2 1'
DATA_FILES='flash.par
amr_runtime_parameters
helm_table.dat
coldwd_mchandra_damped.dat
nse_dens_emq_table.txt
nse_pres_hmq_table.txt
SpeciesList.txt'

echo "Starting Cobalt job script"
for i in ${THREAD_EXPERIMENTS}; do
  echo "Running job with ${i} OpenMP threads"
  echo ""
  echo ""
  echo ""

  mkdir ${i}
  cd ${i}
  ln -s ../${FLASH_BINARY}
  for j in ${DATA_FILES}; do
      ln -s ../${j}
  done

  cobalt-mpirun -nofree -mode smp -np ${NODES} \
      -env "XLSMPOPTS=stack=16000000 OMP_NUM_THREADS=${i}" `pwd`/${FLASH_BINARY}
  echo "Return code is " $?
  cd ..
done

echo ""
echo ""
echo "Freeing partition"
cobalt-mpirun -free wait

#bg-listblocks
#If any blocks still allocated:
#mpirun -partition ANL-R00-M1-512 -free wait
