#!/usr/bin/gnuplot
#
# Showcasing the viridis colormap
#
# AUTHOR: Hagen Wierstorf
# gnuplot 5.0 patchlevel 3

reset

# wxt
#set terminal wxt size 350,262 enhanced font 'Verdana,10' 
# png
set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
set output 'virus.png'

unset key

# border
set style line 11 lc rgb '#808080' lt 1
set border 3 front ls 11
set tics nomirror out scale 0.75

# colorbar
# disable colorbar tics
set cbtics scale 0



plot 'initial_lattice.dat'  matrix with image

