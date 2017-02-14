#!/bin/bash
#
# This script helps us debug problems on Blue Gene machines.  It parses the
# lightweight core files to find the instruction which caused the unhandled
# signal and then converts this instruction into a source file and line number.
#
# Example usage
#
# The .error file contains the MPI rank which had the unhandled signal.
#  EAS-20040-31371-128:8963:ibm.runjob.client.Job: terminated by signal 11
#  EAS-20040-31371-128:8963:ibm.runjob.client.Job: abnormal termination by signal 11 from rank 240
#
# We would run this script as follows
# > ./core_signal_instruction.sh core.240
# (The script assumes that the flash4 binary is in the current working directory)
#
# core.240
# ***FAULT Encountered unhandled signal 0x0000000b (11) (SIGSEGV)
# /bgsys/drivers/DRV2012_0131_2135/ppc64-rhel60/toolchain/gnu/glibc-2.12.2/gmon/../sysdeps/posix/profil.c:82
# source code at line 82:
#     ++samples[i];
#
# You may wonder why I have this script when there is bgq_stack
# > /soft/debuggers/scripts/bin/bgq_stack ./flash4 core.240 
# In the example above bgq_stack shows the application stack traces and does not
# show profil.c in which we encountered the SIGSEGV.  This is because profil.c
# is executed on an interupt and is not part of the application.
#

BINARY='./flash4'

if [ $# -lt 1 ]; then
    echo "Usage $0 core.*"
    exit 1
fi

if [ ! -e ${BINARY} ]; then
    echo "The binary ${BINARY} does not exist"
    exit 2
fi

for i in "$@"; do
    signal=$(grep 'FAULT Encountered unhandled signal' ${i})
    if [ $? -eq 0 ]; then
	echo -e "${i}\n${signal}"
	inst=$(awk -F 0x '/While executing instruction/ { print $2 }' ${i})
	file_and_line=$(addr2line ${inst} -e ${BINARY})
	echo "${file_and_line}"

	if [ "$file_and_line" != "??:0" ]; then
	    file=$(echo ${file_and_line} | awk -F : '{print $1}')
	    if [ -f "${file}" ]; then
		line=$(echo ${file_and_line} | awk -F : '{print $2}')
		source=$(sed -n ${line}p ${file})
		echo -e "source code at line ${line}:\n${source}"
	    fi
	fi
	echo
    fi
done
