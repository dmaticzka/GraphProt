#!/bin/bash

LC_PERF=$1
PLOT_PS=$2

gnuplot -p << EOF
ftr(x)=atr-btr/(x+ctr)
fts(x)=ats-bts/(x+cts)
fit ftr(x) '$LC_PERF' u 1:2  via atr,btr,ctr
fit fts(x) '$LC_PERF' u 1:3  via ats,bts,cts
plot '$LC_PERF' u 1:2 t "" w p lt 1, '' u 1:3 t "" w p lt 2, ftr(x) t "Training measure" w l lt 1 lw 3, fts(x) t "Test measure" w l lt 2 lw 3
set terminal svg
set out '$PLOT_PS'
replot
EOF

