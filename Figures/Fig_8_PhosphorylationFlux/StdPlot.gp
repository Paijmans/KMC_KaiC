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
set yrange [0:0.5]


#filename u 1:5 w l ls 6 t 'ACII',\
#filename u 1:(($7+$8)/12) w l ls 7 t 'ATP frac',\
#

set key 

p filename u 1:(1-$22) w l ls 4 t 'inactive',\
filename u 1:($12/6) w l t 'pT' lt 2 lw 3 lc rgb 'green',\
filename u 1:($14/6) w l t 'pS' lt 2 lw 3 lc rgb 'red',\
filename u 1:($13/6) w l t 'pD' lt 2 lw 3 lc rgb 'blue',\
filename u 1:($7+$8)/12 w l ls 9 t 'ACI',\
filename u 1:($5) w l ls 6 t 'ACII',\
filename u 1:2 w l t 'p' lw 3 lc rgb 'orange',\
filename u 1:3 w l t 'Afree' lw 3 



#p filename u 1:20 w l ls 4 t 'active',\
#filename u 1:($12/6) w l t 'pT' lt 2 lw 3 lc rgb 'green',\
#filename u 1:($14/6) w l t 'pS' lt 2 lw 3 lc rgb 'red',\
#filename u 1:($13/6) w l t 'pD' lt 2 lw 3 lc rgb 'blue',\
#filename u 1:($7/6) w l ls 9 t 'CIATP',\
#filename u 1:4 w l ls 6 t 'ACI'

#filename u 1:((6-$8)/6) w l ls 7 t 'CIIADP'
#filename u 1:3 w l ls 9 t 'Afree',\
#filename u 1:($11/6) w l t 'pU' lt 2 lw 3 lc rgb 'pink',\
#filename u 1:(($7+$8)/12) w l ls 9 t 'ATPtot'

