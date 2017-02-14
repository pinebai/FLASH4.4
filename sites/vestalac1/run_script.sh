#!/bin/bash -x
qsub -A TurbNuclComb_esp -n 32 --mode script -t 1:00:00 ./cobalt_script.sh
