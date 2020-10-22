set style line 1 linewidth 3
set style line 2 linewidth 3
set style line 3 linewidth 3
set style line 4 linewidth 3
set style line 5 linewidth 3
set style line 6 linewidth 3
set style line 7 linewidth 3 lt 3 lc rgb 'black'
set style line 8 linewidth 3 lt 3 lc rgb 'purple'
set style line 9 linewidth 3 lt 3 lc rgb 'gray'

filename='flux.dat'

set key 

p filename u 1:6 w l t 'UT' lt 2 lw 3 lc rgb 'black',\
filename u 1:7 w l t 'TD' lt 2 lw 3 lc rgb 'green',\
filename u 1:8 w l t 'SD' lt 2 lw 3 lc rgb 'blue',\
filename u 1:9 w l t 'US' lt 2 lw 3 lc rgb 'red',\
filename u 1:($3+$4+$5)/10 w l,\
0

#p filename u ($6+$7):($8+$9) t 'UT' lt 2 lw 3 lc rgb 'black'
