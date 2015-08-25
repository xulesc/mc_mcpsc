# plotroc.p

## Note the solid makes lines solid
set term postscript eps enhanced color solid

set out "sc_psc_roc.eps"

set key bottom right

set key font ",16"
set xtics font ",14"
set xlabel "False Positive Rate (FPR)" font ",18"
set ylabel "True Positive Rate (TPR)" font ",18"

plot 'perf/mcpsc_perf.dat' using 9:8 title 'MCPSC' lc rgb "black" with lines, \
 'perf/tmalign_perf.dat' using 9:8 title 'TMALIGN' lc rgb "red" with lines, \
 'perf/ce_perf.dat' using 9:8 title 'CE' lc rgb "green" with lines, \
 'perf/gralign_perf.dat' using 9:8 title 'GR-ALIGN' lc rgb "blue" with lines, \
 'perf/fast_perf.dat' using 9:8 title 'FAST' lc rgb "cyan" with lines, \
 'perf/usm_perf.dat' using 9:8 title 'USM' lc rgb "purple" with lines;
 

set out "sc_psc_pr.eps"

plot 'perf/ce_perf.dat' using 7:6 title 'CE' with lines, \
 'perf/fast_perf.dat' using 7:6 title 'FAST' with lines, \
 'perf/gralign_perf.dat' using 7:6 title 'GR-ALIGN' with lines, \
 'perf/tmalign_perf.dat' using 7:6 title 'TMALIGN' with lines, \
 'perf/usm_perf.dat' using 7:6 title 'USM' with lines, \
 'perf/mcpsc_perf.dat' using 7:6 title 'MCPSC' with lines;



 
 