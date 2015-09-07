# plot.p

set term postscript eps enhanced color solid

set out "sc_psc_bar_cath.eps"

set style data histogram
set style histogram rowstack gap 1
set style fill solid border -1
set boxwidth 0.9
set key autotitle columnheader
set key outside right center vertical
set xtics nomirror rotate by -30 font ",12"
set ytics font ",20"
set bmargin 6
set key font ",24"

C = "#99ffff"; Cpp = "#4671d5"; Java = "#ff0000"; Python = "#f36e00"

set style histogram rowstacked title offset 0,-1 

plot newhistogram "{/*1.5 Class 1}" lt 1, 'data.dat' using 3:xtic(1) lc rgb C, '' u 4 lc rgb Cpp, \
 newhistogram "{/*1.5 Class 2}" lt 3, '' u 5:xtic(1) lc rgb Python , '' u 7 lc rgb Cpp, \
 newhistogram "{/*1.5 Class 3}" lt 5, '' u 8:xtic(1) lc rgb  Python, '' u 9 lc rgb C;
