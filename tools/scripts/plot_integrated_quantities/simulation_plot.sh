#!/usr/bin/env bash
#
# This script uses gnuplot to create integrated quantity plots for 1
# to N data files.  The output is saved in a single postscript file
# which is then converted to pdf if the ps2pdf command (belongs to the
# ghostscript package) is in $PATH.
#
#
# The pdf file is generally more useful than the postscript file, but
# I keep the postscript file around because there are rarely pdf
# viewers installed on supercomputers.  The postscript file can be
# viewed with gs (belongs to the ghostscript package) which is
# installed on most machines including intrepid and mira.
#
#
# The input parameters are obtained by sourcing a script in the
# current working directory named plot_parameters.sh.  The parameters
# are
#  FILE1: 1st .dat file (compulsory)
#  FILE2: 2nd .dat file (optional)
#  FILE3: 3rd .dat file (optional)
#  FILE1_LABEL: 1st .dat file label appearing in the legend (optional)
#  FILE2_LABEL: 2nd .dat file label appearing in the legend (optional)
#  FILE3_LABEL: 3rd .dat file label appearing in the legend (optional)
#  PLOT_NAME: Name of the output plot file - no extension needed (optional)
#  PLOT_TITLE: Title appearing at the top of each plot (optional)
#
# The plot_parameters.sh file is not mandatory; plots can be created
# by directly setting the bash variables, e.g.
#  $ FILE1=bgp_rtflame.dat simulation_plot.sh
#  $ FILE1=bgp_rtflame.dat FILE2=bgq_rtflame.dat simulation_plot.sh
#
#
# We process the data in the input files before we plot the data with
# gnuplot (the input files are untouched).  The processing is as follows
#   1. Remove any binary characters from the data files.  We found
#      binary characters in a BG/Q data file after an I/O node panic.
#      [tr command].
#   2. Remove the first set of a data points whenever there is a
#      restart.  This is because some quantities such as "Burning
#      rate" show 0.0 after a restart because no previous value is
#      available.  [sed command].
#   3. Sort the data.  This is because the input file can be a
#      concatenation of multiple integrated quantities files.  In the
#      event of a run failure there can be values in a specific
#      integrated quantities file at a later time than the last
#      checkpoint time.  If we restarted from this checkpoint file
#      there will be duplicate data points and so sorting is a quick
#      and dirty way to prevent gnuplot joining data points with lines
#      in a strange way [sort command].
#
#
# This script is known to work with gnuplot >= 4.0.

set -e

PLOT_PARAMETERS="plot_parameters.sh"
ERR_NO_GNUPLOT=-1
ERR_NO_FILE=-2
ERR_BAD_FILE=-3
ERR_UNDEFINED_VAR=-4

if ! command -v gnuplot >/dev/null 2>&1; then
    echo "No gnuplot in PATH... exiting"
    exit ${ERR_NO_GNUPLOT}
fi

if [ -f ${PLOT_PARAMETERS} ]; then
    echo "Found ${PLOT_PARAMETERS} file... sourcing file"
    source ${PLOT_PARAMETERS}
fi
[ -z "${PLOT_NAME}" ] && PLOT_NAME="integrated_quantities"
[ -z "${PLOT_TITLE}" ] && PLOT_TITLE="FLASH integrated quantities as a function of time"


i=1
while :; do
    # Use an explanation mark inside the braces to retrieve the value:
    # http://www.gnu.org/software/bash/manual/bashref.html#Shell-Parameter-Expansion
    file_variable="FILE${i}"
    file="${!file_variable}"

    if [ -z "${file}" ]; then
	break
    else
	echo "${file_variable} is ${file}"

        # Check that the data file exists.
	if [ ! -f "${file}" ]; then
	    echo "${file} does not exist!"
	    exit ${ERR_NO_FILE}
	fi

	# Check that the data file has at least 2 columns and that it
	# has the same number of columns as all other data files.
	ncol=$(awk '{ if (NR==2) {print NF; exit} }' "${file}")
	if [ $ncol -lt 2 ]; then
	    echo "There must be >= 2 columns in ${file}!)"
	    exit ${ERR_BAD_FILE}
	fi
	if [ "${i}" == "1" ]; then
	    ncol_file1=${ncol}
	else
	    if [ $ncol -ne $ncol_file1 ]; then
		echo "Different number of columns in FILE1 and ${file}!"
		exit ${ERR_BAD_FILE}
	    fi
	fi

	files[i]="${file}"
	let i=i+1
    fi
done

# Add || true because i-1 may evaluate to 0:
# http://mywiki.wooledge.org/BashFAQ/105/Answers
let numfiles=i-1 || true
if [ $numfiles -lt 1 ]; then
    echo "You must specify at least one data file (define FILE1)!"
    exit ${ERR_UNDEFINED_VAR}
fi



# Use the cut command to extract the x and y plot labels.
# We assume the labels are the same in all data files.
FILE_COLUMNS[0]="time (seconds)" # Looks better than '#time'
let nquant=ncol-1
offset=2; width=26 # Fortran format string is format (2x,50(a25, :, 1X))
col=1
for i in $(seq 1 ${nquant}); do
    let start=(i*width)+offset+1
    let end=start+width-1
    ylabel=$(cut -c ${start}-${end} "${FILE1}" | head -1 | sed 's/ *$//g')
    if [ -n "${ylabel}" ]; then
	FILE_COLUMNS[${col}]="${ylabel}"
	echo "Integrated quantity ${col} has title ${FILE_COLUMNS[${col}]}"
	let col+=1
    else
	echo "Integrated quantity ${col} has blank title... break from loop"
	break
    fi
done

if [ ${#FILE_COLUMNS[*]} -lt 2 ]; then
    echo "We need at least two columns to create a 2D plot... exiting"
    exit ${ERR_BAD_FILE}
fi



# Remove binary characters from the data files and remove the first set of
# data points after a FLASH restart.  Also find a time range that is common
# across all data files.
for i in $(seq 1 ${numfiles}); do
    file="${files[i]}"
    tmpfile="$(mktemp --tmpdir=.)"

    tr -cd '\11\12\15\40-\176' < "${file}" | \
	sed '/# simulation restarted/,+1d' | sort -g > "${tmpfile}"
    tmin_file=$(awk '{ if (substr($0,3,1) != "#") {print $1; exit} }' "${tmpfile}")
    tmax_file=$(awk 'END {print $1; exit}' "${tmpfile}")
    echo "FILE${i} has data from ${tmin_file} to ${tmax_file} seconds"

    if [ "${i}" == "1" ]; then
	tmin=${tmin_file}
	tmax=${tmax_file}
    else
	tmin=$(echo "${tmin} ${tmin_file}" | \
	    awk '{ p=$1; q=$2; if (p>q) {tmin=p} else {tmin=q}; print tmin}')
	tmax=$(echo "${tmax} ${tmax_file}" | \
	    awk '{ p=$1; q=$2; if (p>q) {tmax=q} else {tmax=p}; print tmax}')
    fi

    tmpfiles[i]="${tmpfile}"
done

echo -e "          plot from ${tmin} to ${tmax} seconds\n"
if [ -z ${tmin} -o -z ${tmax} ]; then
    echo "Invalid time range!  Exiting"
    for file in "${tmpfiles[@]}"; do rm ${file}; done
    exit ${ERR_UNDEFINED_VAR}
fi



# Generate the gnuplot plotting script.
GNUPLOT_SCRIPT="$(mktemp --tmpdir=.)"
cat > "${GNUPLOT_SCRIPT}" << EOF
#!/usr/bin/env gnuplot
set term postscript noenhanced color
set title "${PLOT_TITLE}"
set output "${PLOT_NAME}.ps"
set key outside top
set rmargin 4
set xlabel "${FILE_COLUMNS[0]}"
set xrange [${tmin}:${tmax}]
EOF

for col in ${!FILE_COLUMNS[*]}; do
    if [ $col -gt 0 -a $col -lt $ncol ]; then
	ylabel_cmd="set ylabel \"${FILE_COLUMNS[${col}]}\""

	i=1
	plot_cmd=""
	let gnuplot_col=col+1
	for file in "${tmpfiles[@]}"; do
	    if [ "${i}" == "1" ]; then
		line_prefix="plot"
	    else
		line_prefix=","
	    fi

	    file_label_variable="FILE${i}_LABEL"
	    file_label="${!file_label_variable}"
	    [ -z "${file_label}" ] && file_label="FILE${i}"

	    plot_cmd+="${line_prefix} \"${file}\" using 1:${gnuplot_col} "\
"title \"${file_label}\" with lines lw 4 lt ${i}"

	    let i=i+1
	done

	# Write the plot and ylabel commands to the script.
	echo "$ylabel_cmd" >> "${GNUPLOT_SCRIPT}"
	echo "$plot_cmd" >> "${GNUPLOT_SCRIPT}"
    fi
done



# Run the script, clean up and convert the postscript file into pdf.
echo "Executing gnuplot script..."

chmod +x "${GNUPLOT_SCRIPT}"
gnuplot "${GNUPLOT_SCRIPT}"
rm "${GNUPLOT_SCRIPT}"
echo "The postscript output is saved in ${PLOT_NAME}.ps"

for file in "${tmpfiles[@]}"; do rm ${file}; done

if command -v ps2pdf >/dev/null 2>&1; then
    echo "Found ps2pdf (belongs to ghostscript package)... converting to pdf"
    ps2pdf "${PLOT_NAME}.ps"
    echo "The pdf output is saved in ${PLOT_NAME}.pdf"
else
    echo "Did not find ps2pdf (belongs to ghostscript package)... continuing"
fi

echo "Done."
exit 0
