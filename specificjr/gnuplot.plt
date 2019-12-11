set terminal pdf
set output "rj.pdf"
set title "Specific heat"
set xlabel "T"
set ylabel "Cp"
plot "op.txt" using 1:2 with lines title "cp"