# plot-save.plt
set terminal pngcairo size 2000, 2000 fontscale 4
set style fill solid border -1
set key off


binwidth=1.5
bin(x,width)=width*floor(x/width)


set title "A * A'  Eigenvalues"
set output "hist-mul.png"
plot 'data-p2-mul.txt' using (bin($1,binwidth)):(1.0) smooth freq with boxes

clear

set title "A + A'  Eigenvalues"
set output "hist-sum.png"
plot 'data-p2-sum.txt' using (bin($1,binwidth)):(1.0) smooth freq with boxes



