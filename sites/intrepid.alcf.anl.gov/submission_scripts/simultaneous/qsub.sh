#!/bin/bash
qsub -q prod -A TurbNuclComb_esp -t 720 -n 8192 --mode script batch_submission_script
