#!/usr/bin/gnuplot -persist

set title "Quality of finite-difference solution when using \n\
16^{3}, 32^{3}, 64^{3} and 128^{3} grid points"

set xlabel "log (delx)" font "arial,20"
set ylabel "log (Linf error)" font "arial,20"
set key 0.002, 0.02
set pointsize 6

set rmargin 4
set logscale x
set logscale y
set format x "10^{%L}"
set format y "10^{%L}"

#We are plotting the Linfinity error as a function of cell spacing.  
#Cell spacing is 1/Gridsize (grid has the same length along each dimension).
#We plot delx^2 to demonstrate that the Linfinity error is decreasing at 
#a rate proportional to delx^2.

plot "pfft_results.txt" using (1/$1):2 title "pfft " with linespoints, \
"pfft_results.txt" using (1/$1):((1/$1)*(1/$1)) title "delx^{2}" with linespoints
