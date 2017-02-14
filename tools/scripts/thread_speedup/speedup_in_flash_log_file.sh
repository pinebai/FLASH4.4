#!/bin/bash -x
#
# This script enters directories ${DIRS} and extracts timing information from
# the flash log file ${LOGFILE}.  The name ${DIR} indicates the number of
# threads used for the run, e.g. 1,2,4.  A single column of timing information
# ${TIMER_COLS} is stored from each log file and written to ${OUTPUT_FILE}.
#
# This script assumes a directory named 1 is present.
#

LOGFILE="flash.log"
#LOGFILE="rtflame_256_s9.log"
OUTPUT_FILE="time_n_threads.txt"

DIRS="1 2 4"
TIMER_STR="accounting unit                  max/proc"
FILE_A="file.tmp.a"
FILE_B="file.tmp.b"
FILE_C="file.tmp.c"
TIMER_ROWS="200"
TIMER_COLS="35-46" #Use this for max/proc
#TIMER_COLS="62-73" #Use this for avg/proc

#Check that the run directories exist.
for i in ${DIRS}; do
    if [ ! -d "$i" ]; then
	echo "Directory $i does not exist!"
	exit -1
    fi
done

#Copy the labels into FILE_A
grep "${TIMER_STR}" -A ${TIMER_ROWS} 1/${LOGFILE} | cut -c 1-30 > ${FILE_A}

for i in ${DIRS}; do
    #Extract the raw timing information from the logfile into FILE_B.
    grep "${TIMER_STR}" -A ${TIMER_ROWS} ${i}/${LOGFILE} | cut -c ${TIMER_COLS} > ${FILE_B}
    paste ${FILE_A} ${FILE_B} > ${FILE_C}

    #Copy the current output file into FILE_A.
    cp ${FILE_C} ${FILE_A}
done

rm ${FILE_A}
rm ${FILE_B}

end_line=$(awk '/===========/{print NR-1}' ${FILE_C})
echo "                                1 thread        2 threads       4 threads    " > ${OUTPUT_FILE}
head "-${end_line}" ${FILE_C} >> ${OUTPUT_FILE}

rm ${FILE_C}
