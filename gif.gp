set term gif animate delay 15

set output 'Gif_F.gif'

set ylabel 'Probabilidad' font 'Arial,15'

set xlabel 'Posición (x)' font 'Arial,15'

unset key

set xra [0:1]
set yra [-0.00275:0.00275]

do for [i=0:150:1]{plot 'prob_gif.dat' index i w l}