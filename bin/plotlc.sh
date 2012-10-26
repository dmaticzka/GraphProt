#!/bin/bash

LC_PERF=$1
PLOT_PS=$2

gnuplot -p << EOF
set terminal png
set grid
set xlabel "Training set size"
set ylabel "Area Under ROC Curve"
set out '$PLOT_PS'
ftr(x)=atr-btr/(x+ctr)
fts(x)=ats-bts/(x+cts)
fit ftr(x) '$LC_PERF' u 1:9  via atr,btr,ctr
fit fts(x) '$LC_PERF' u 1:17  via ats,bts,cts
plot '$LC_PERF' u 1:9 t "" w p lt 1, '' u 1:17 t "" w p lt 2, ftr(x) t "Train" w l lt 1 lw 2, fts(x) t "Test"  w l lt 2 lw 2
EOF

rm fit.log
exit 0;

