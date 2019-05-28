#!/bin/bash
#this script simply plot the two curves
#ponychen 20190528

gnuplot -persist << EOF
set xlabel "engineering strain" 
set ylabel "engineering stress/GPa"
set title "engineering stress strain curve"
plot "engineering_stress_strain.all"  w l lw 6 lc 8
EOF

gnuplot -persist << EOF
set xlabel "true strain"
set ylabel "true stress/GPa"
set title "true stress strain curve"
plot "true_stress_strain.all"  w l lw 6 lc 9
EOF
