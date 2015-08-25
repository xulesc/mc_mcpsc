# plot.p

set term postscript eps enhanced color solid

set out "sc_psc_bar.eps"

set style data histogram
set style histogram rowstack gap 1
set style fill solid border -1
set boxwidth 0.9
set key autotitle columnheader
set key outside right center vertical
set xtics nomirror rotate by -30 font ",12"
set bmargin 6

set style histogram rowstacked title offset 0,-1 

plot newhistogram "Class 1" lt 1, 'data.dat' using 3:xtic(1) , '' u 4 , \
 newhistogram "Class 2" lt 3, '' u 5:xtic(1) , '' u 7 , \
 newhistogram "Class 3" lt 5, '' u 8:xtic(1) , '' u 9 ;
