set term png

set output 'Prob_F.png'

set ylabel 'Probabilidad' font 'Arial,15'

set xlabel 'Posición (x)' font 'Arial,15'

set key box
set key left

set xra [0:1]
set yra [-0.00275:0.00275]

plot 'prob_png.dat' index 0 w l lw 2 t ' t=0 seg ', 'prob_png.dat' index 1 w l lw 2 t ' t=1E-3 seg ', 'prob_png.dat' index 2 w l lw 2 t ' t=2E-3 seg ', 'pot.dat' w l lt 0 lw 2 t ' Potencial '