set term qt

plot "OUTPUT_FILES/AA.S0001.BXY.semd" t 'Ux Heterogeneous' w l lc 1

set xlabel "Time (s)"
set ylabel "Amplitude of displacement component (m)"

set xrange [-0.15:30]
set yrange [-0.06:0.06]

pause -1 "Hit any key..."


