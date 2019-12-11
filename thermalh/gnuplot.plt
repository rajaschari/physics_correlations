set terminal pdf
set output "thermal.pdf"
set title "Thermal Hall coefficient"
set xlabel "T"
set ylabel "TH"
plot "op.txt" using 1:2 with lines title ""