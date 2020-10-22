set terminal epscairo

datafile="Statespace_flux.hist"

set output 'Statespace_flux.eps'

set cbtics out
unset key
set notitle

#set palette rgb 33,13,10;
#set palette rgb 7,5,15

set palette defined ( 0 "#7B0099",\
                      1 "#000090",\
                      2 "#000fff",\
                      3 "#0090ff",\
                      4 "#0fffee",\
                      5 "#90ff70",\
                      6 "#ffee00",\
                      7 "#ff7000",\
                      8 "#ee0000",\
                      9 "#7f0000")




set xlabel '#pT'
set ylabel '#pS'

#set xrange  [  -6.5 : 6.5 ] noreverse nowriteback\
set xrange  [  -0.5 : 6.5 ] noreverse nowriteback
set yrange  [  -0.5 : 6.5 ] noreverse nowriteback
set cbrange [  0 : 0.15 ] noreverse nowriteback                                       
set cblabel "#Hex(#pT,#pS)" offset 1,0

#plot datafile using 1:2:3 with image,\
#'' every 1 u 1:2:(0.1*$4/sqrt($4*$4+$5*$5)):($3>0.001?0.1*$5/sqrt($4*$4+$5*$5):1/0) with vectors arrowstyle 31

plot datafile using 1:2:3 with image,\
'' every 1 u 1:2:(100*$4):($3>0.001?100*$5:1/0) with vectors filled lc rgb 'white'

#plot datafile using 1:2:3 with image,\
#'' every 1 u 1:2:($4):($3>0.001?$5:1/0) with vectors arrowstyle 31
