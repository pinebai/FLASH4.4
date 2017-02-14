#!/bin/bash
#
# This script compares Chombo HDF5 checkpoint files.  It is used because sfocu
# does not understand Chombo file layout yet.  It compares the checkpoint
# files using h5diff and writes "SUCCESS" or "FAILURE" to stdout depending
# on the h5diff exit status code.  The "SUCCESS" text is required by the
# FLASH test suite program (flashTest.py).
#
# Optional arguments:
#   -t: The full filename of the comparison tool (h5diff).
#
# The test suite will use this script in place of sfocu if the config file
# is modified as follows:
#  pathToSfocu: /full/path/to/cmp_chombo.sh -t /full/path/to/h5diff
#  sfocuScript: <pathToSfocu>
# (The -t optional argument is used in case h5diff is not in $PATH.)
#
# The article "Parsing arguments for your shell script" at
# http://www.linux.com/archive/feed/118031 explains how to handle
# optional arguments.
#

EXPECTED_ARGS=2
BAD_ARG=65

function exit_with_help {
    printf "Usage: ./%s [-t comparison_tool] file1 file2\n" $(basename $0) >&2
    exit ${BAD_ARG}
}


#Assume h5diff is in $PATH by default.  The user may add the -t optional
#argument to specify the full filename of the h5diff executable.
cmp_tool="h5diff"
while getopts 't:' OPTION
do
    case $OPTION in
	t)cmp_tool="$OPTARG"
	    ;;
	?)exit_with_help
	    ;;
    esac
done
shift $(($OPTIND - 1))


#Check that we have two arguments remaining.
if [ $# -ne ${EXPECTED_ARGS} ]
then
    exit_with_help
else
    file1="$1"
    file2="$2"
    ${cmp_tool} -c ${file1} ${file2}
    rc=$?
    if [ ${rc} -eq 0 ] ; then
	echo "SUCCESS"
    else
	echo "FAILURE - h5diff return code is ${rc}"
    fi
    exit ${rc}
fi
