#!/bin/bash
#This script produces a png image that compares jeans analytical 
#solution to the simulation solution over the first 2.5 seconds.
#Depends on gfortran, gnuplot, sed, awk, head.

DATA_FILE="jeans.dat"
SIM_PLOT_FILE="JeansSimulation.dat"
ANA_PLOT_FILE="JeansAnalytical.dat"
GNUPLOT_SCRIPT="Compare.plt"

#Find the internal energy at time t=0.
u0_=$(sed '1,2d' ${DATA_FILE} | awk '{print $6'} | head -1)

#Print the simulation solution.
sed '1,2d' ${DATA_FILE} | awk -v "u0=${u0_}" '{ ut=$8; ur=(ut-u0); printf "%.16e, %.16e, %.16e, %.16e\n", $1,$7,ur,$9*10.0 }' > $SIM_PLOT_FILE

#Print the analytical solution.
gfortran JeansAnalytical.F90 -o JeansAnalytical
./JeansAnalytical

cat > "$GNUPLOT_SCRIPT" << EOF
#!/usr/bin/gnuplot
set terminal png enhanced size 1024,768
set output "jeans_validate.png"
set key outside
set xrange [0:2.5]
set xlabel "Time (sec)"
set ylabel "Energy (erg)"
plot "${ANA_PLOT_FILE}" using 1:2 title "Kinetic Energy (analytical) T(t)" with lines lt 3, \
"${ANA_PLOT_FILE}" using 1:3 title "Internal Energy (analytical) U(t)-U(0)" with lines lt 4, \
"${ANA_PLOT_FILE}" using 1:4 title "Potential Energy (analytical) W(t)*10" with lines lt 1, \
"${SIM_PLOT_FILE}" using 1:2 title "Kinetic Energy (simulation) T(t)" with points lt 3 pt 1 ps 2.5, \
"${SIM_PLOT_FILE}" using 1:3 title "Internal Energy (simulation) U(t)-U(0)" with points lt 4 pt 2 ps 2.5, \
"${SIM_PLOT_FILE}" using 1:4 title "Potential Energy (simulation) W(t)*10" with points lt 1 pt 3 ps 2.5
#pause -1
EOF

chmod +x "${GNUPLOT_SCRIPT}"
./"${GNUPLOT_SCRIPT}"
rm -f "${SIM_PLOT_FILE}" "${ANA_PLOT_FILE}" "${GNUPLOT_SCRIPT}" JeansAnalytical JeansAnalytical.o
