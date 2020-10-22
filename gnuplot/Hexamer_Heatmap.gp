datafile="Hexamers.dat"

set cbtics out
unset key
set notitle

set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")




set xlabel 'Time'
set ylabel 'Hexamer #'

set xrange  [  0 : 720 ] noreverse nowriteback
set yrange  [  0 : 120 ] noreverse nowriteback
#set cbrange [  0 : 6 ] noreverse nowriteback                                       
set cblabel "Hexamer state" offset 1,0

plot datafile using 1:2:3 with image
