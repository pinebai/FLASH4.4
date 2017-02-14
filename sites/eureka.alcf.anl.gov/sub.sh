#!/bin/bash
#Pass 'mx' or 'nomx' to run_script.sh
qsub -A SupernovaVandV -n 1 -t 5 --mode script ./run_script.sh 'mx'
