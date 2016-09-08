set term pdfcairo
set output "1D_HO_propagator.pdf"
set size 1, 1
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "1D Harmonic Oscillator Propagator"
set ylabel "Propogator Amplitude"
set xlabel "x" 
plot "1D_HO_propagator.txt" using 1:2 title "path integral calculation" with lines, "1D_HO_propagator.txt" using 1:3 title "exact" with lines, 
set output
