# plot-save.plt
set terminal pngcairo size 8000, 8000 fontscale 10 pointscale 1
set size square
set output "sin-plot.png"
set xrange[0:5]
set yrange[-1.5:1.5]
set title "y = Sin(2x)"

plot "data-p1-sin.txt" with linespoints pt 30 pointsize 10 title ""

clear

set xrange[-2.5:2.5]
set yrange[-2.5:2.5]
set output "surface-plot.png"
set title "z = Tanh(x*y)"
set hidden3d

splot "data-p1-surface.txt" with lines lc rgb '#b90046' lw 1.5 title ""



