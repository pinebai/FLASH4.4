#!/bin/bash
#
# This script uses HPCToolkit differential profiling to find
# scalibility bottlenecks.  It does this by combining multiple
# execution profiles from experiments run on different numbers
# of processors into a single profile file.  It takes each
# processor value P in PROCS array, locates flash.par.P and then
# runs flash4 on P processors.  It will perform all the steps to
# create a single HPCToolkit profile named experiment.xml which
# can be viewed using hpcviewer.
#
# This script should be added to your $PATH.  The script assumes
# that hpcrun and mpiexec are in $PATH.
#
# The script accepts the following command line options:
# -n : number of MPI processes (PROCS)
#      * for a single experiment use e.g. -n 8
#      * for multiple experiments use e.g. -n '1 4' (note the quotes)
#      (defaults to '1 4', i.e. 1 and 4 MPI processes)
# -s : directory containing FLASH source (OBJDIR)
#      (defaults to current directory)
# -r : directory where flash4, flash.pars and other run files exist (RUNDIR)
#      (defaults to current directory)
#
#
# Example 1 - Obtain a profile for a FLASH run with 1 MPI process.
#             We launch the script from FLASH objdir directory.
#             (it will use flash.par.1 parfile)
#
#   differential_profile.sh -n 1
#
#
# Example 2 - Obtain profiles for FLASH runs with 2 and 4 MPI processes.
#             We launch the script from a run directory so need to specify
#             the path to FLASH objdir directory.
#             (it will use flash.par.2 and flash.par.4 parfiles)
#
#   differential_profile.sh -n '2 4' -s /path/to/flash/objdir
#
#
# Example 3 - Obtain a profile as in Example 1 but on Eureka and
#             with 64 MPI processes.
#             (it will use flash.par.64 parfile)
#
#   qsub -A SupernovaModels -n 8 -t 180 --mode script ./script.sh
#
# where script.sh is:
#   mpdboot -n 8 -f $COBALT_NODEFILE
#   export PATH=$PATH:/home/cdaley/eureka/software/hpctoolkit/5.1.0/bin
#   differential_profile.sh -n 64
#

PROCS='1 4'
OBJDIR=$(pwd)
RUNDIR=${OBJDIR}

BINARY='flash4'
ALL_RUNFILES='amr_runtime_parameters SpeciesList.txt eos_helm.dat'
PROFILE_DIR='differential_profile'
FILE=${PROFILE_DIR}/'experiment.xml'
FILE_BAK=${PROFILE_DIR}/'experiment.xml.bak'

# Command line options for additional flexibility.
while getopts "n:r:s:" opt
do
  case $opt in
      n) PROCS=$OPTARG;;
      s) OBJDIR=$OPTARG;;
      r) RUNDIR=$OPTARG;;
  esac
done


#Enter the run directory and compose a list of all the runtime files.
cd ${RUNDIR}

runfiles=''
for runfile in ${ALL_RUNFILES}; do
    if [ -f ${runfile} ]; then
	runfiles=${runfiles}' '${runfile}
    fi
done


#Run experiments with i MPI processes.
for i in ${PROCS}; do
    parfile="flash.par.${i}"
    if [ -f ${parfile} ]; then
	echo "Found parfile ${parfile} for ${i} processor run"
	mkdir ${i}
        cd ${i}
	ln -s ../${BINARY}
	ln -s ../${parfile}
	for j in ${runfiles}; do
	    ln -s ../${j}
	done
	hpcrun -o profile mpiexec -n ${i} ./${BINARY} -par_file ${parfile}
	cd ..
    else
	echo "Did not find parfile ${parfile} for ${i} processor run"
    fi
done


#Create a profile by combining all raw measurements from directories i.
s=''
for i in ${PROCS}; do
    if [ -d "${i}/profile" ]; then
        #The 6 zeros represent MPI rank 0.  On Eureka all profile
        #files have MPI rank 0 for some reason.  For now we take the
        #first profile file found - this is not a good solution.
        #The next 3 zeros represent the thread ID - even when running
        #single threaded FLASH we may get several thread IDs
        #because some MPI implementations use helper threads.
	f=$(ls -1 ${i}/profile/${BINARY}-000000-000-*.hpcrun | tail -1)
	s=${s}' '${f}
    fi
done

if [ ! -z "${s}" ]; then
    echo "Running hpcstruct on ${BINARY} binary"
    hpcstruct ${BINARY}

    echo "About to run hpcprof-mpi with measurement files ${s}"
    proc=${PROCS[@]: -1}
    mpiexec -n ${proc} hpcprof-mpi --name ${BINARY}_scaling -o ${PROFILE_DIR} -S ${BINARY}.hpcstruct -I${OBJDIR} ${s}


    #Now that we have the experiment.xml profile we edit the XML
    #to reflect the number of processors used.  We do this because
    #the measurements shown in hpcviewer are labelled 1,2,3 which is
    #not very clear.
    #See https://mailman.rice.edu/pipermail/hpctoolkit-forum/2011-May.txt.gz
    #which has inspired me to start hacking at the raw XML.
    echo "Editing the raw XML"
    cp ${FILE} ${FILE_BAK}
    database=1
    for i in ${PROCS}; do
	if [ -d "${i}/profile" ]; then
	    sed -i "/<Metric i=/ s/n=\"${database}./n=\"${i}PE./" ${FILE}
	    let database=database+1
	fi
    done

    #I find that hpcviewer displays too many metrics by default.
    #Edit the XML to turn off all metrics and then re-add mean
    #inclusive and mean exclusive time.  All metrics are still
    #available and can be turned on within hpcviewer.
    sed -i "/<Metric i=/ s/show=\"1\"/show=\"0\"/" ${FILE}
    sed -i "/Mean ([IE])/ s/show=\"0\"/show=\"1\"/" ${FILE}
else
    echo "No ${BINARY} experiments available"
fi
