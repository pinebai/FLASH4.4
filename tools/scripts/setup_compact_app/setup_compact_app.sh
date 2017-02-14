#!/usr/bin/env bash
#
# This script creates a compact FLASH application.
# Use this script exactly like the standard setup script, e.g.
#
# ./setup_compact_app.sh unitTest/PFFT_PoissonFD -auto -3d 
# -maxblocks=800 +pm4dev -noc -debug -noclobber -parfile=test_paramesh_3d_64cube.par
#
# The additional requirement is a file named UNITS_TO_DELETE_FOR_COMPACT_APP
# in the simulation directory which contains all the units that should be
# removed from our compact application.
#
# The script uses the undocumented --kill-unit setup script option to
# forcefully remove units without caring about REQUESTS or REQUIRES.
# 
# It is then the task of another script named remove_unit_refs.py
# to actually remove all references to the removed units.  The 
# references are the use lines and the subroutine calls.
#
# If the reference removal happens as anticipated you should now be able
# to build the FLASH application as normal.  The object directory in the
# above example is named PFFT_PoissonFD_compact.
#
file="UNITS_TO_DELETE_FOR_COMPACT_APP"
killoption="--kill-unit="
withoutoption="--without-unit="

if [ $# -lt 1 ]; then
    echo "Usage: $(basename $0) <problem-name> [options] [VAR=VALUE]..."
    exit -1
fi

sim="source/Simulation/SimulationMain"
filepath="${sim}/$1"

#obj="$(basename $1)_compact"
obj="$(basename $1)"
partheader="source/Particles/Particles.h"


#Go back to the base level directory
cd ../../..


#If no file found at source/Simulation/SimulationMain/unitTest/Eos/Helmholtz
#then we check source/Simulation/SimulationMain/unitTest/Eos and so on
#until we reach source/Simulation/SimulationMain
while [ ! -f "${filepath}/${file}" ]; do
    echo "${filepath} does not contain ${file}"
    if [[ ${filepath} == ${sim}/* ]]; then
	filepath=$(dirname "${filepath}")
    else
	echo "Cannot find inherited ${file} in ${sim}"
	exit -2
    fi
done
echo "About to remove the units listed in ${filepath}/${file}"


#Construct the setup script options to remove units from our compact app.
removeunits=""
while read line
do
    if [[ ${line} == *flashUtilities* ]]; then
	# This is a special case for flashUtilities which is treated
	# specially by the FLASH setup script.
	removeunits="${removeunits}${withoutoption}${line} "
    else
	removeunits="${removeunits}${killoption}${line} "
    fi
done < "${filepath}/${file}"


#Now call the standard FLASH setup script
setupline="$@ ${removeunits} -objdir=${obj}"
echo -e "The compact app setup line is:\n./setup ${setupline}"
./setup ${setupline}


if [ "$?" -eq "0" ]; then
    #Add a hacky special case for Particles.h header file which is
    #needed by Grid.
    if [[ ${removeunits} == *[Pp]articles* ]]; then
	if [[ ${setupline} == *portable* ]]; then
	    cp $(pwd)/${partheader} $(pwd)/${obj}
	else
	    ln -s $(pwd)/${partheader} $(pwd)/${obj}
	fi
    fi

    #Now call another script to remove references to all removed units
    ./tools/scripts/setup_compact_app/remove_unit_refs.py \
	"${filepath}/${file}" "${obj}"
fi
