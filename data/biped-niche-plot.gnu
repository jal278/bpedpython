set terminal png
set output "biped-niches.png"
set style data lines

set xlabel "Evaluations"
set ylabel "Occupied Niches"
set key bottom right
rows = "3 6"
plot for [i in rows] "biped-niches.dat" u 2:i w l title columnheader(i)

