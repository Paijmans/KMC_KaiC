set style line 1 linewidth 3
set style line 2 linewidth 3
set style line 3 linewidth 3
set style line 4 linewidth 3
set style line 5 linewidth 3
set style line 6 linewidth 3
set style line 7 linewidth 3 lt 3 lc rgb 'black'
set style line 8 linewidth 3 lt 3 lc rgb 'purple'
set style line 9 linewidth 3 lt 3 lc rgb 'gray'

filename='test1.dat'

p filename u 1:($7 + $8)/12 w l ls 1 t 'c_{ATP}',\
filename u 1:($7)/6 w l ls 2 t 'CIATP',\
filename u 1:($8)/6 w l ls 3 t 'CIIATP',\
filename u 1:5 w l t 'ACII' lw 3 lc rgb 'orange',\
filename u 1:($14/6) w l t 'pS' lw 3 lc rgb 'purple',\
filename u 1:20 w l ls 4 t 'active'
#filename u 1:($16/24) w l t 'dATP' ls 5 ,\

