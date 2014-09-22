set terminal png
set output "biped-fitness.png"
set style data lines

set xlabel "Evaluations"
set ylabel "Champion Fitness"
set key bottom right
plot for [i=3:6] "biped-fitness.dat" u 2:i w l title columnheader(i)

