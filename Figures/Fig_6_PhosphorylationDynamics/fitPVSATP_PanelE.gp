set style line 1 linewidth 3
set style line 2 linewidth 3
set style line 3 linewidth 3
set style line 4 linewidth 3
set style line 5 linewidth 3
set style line 6 linewidth 3
set style line 7 linewidth 3 
set style line 8 lw 1 lt 2 lc rgb 'black'

COL1=2
COL2='($7+$8)/12'

unset key

kf=0.50
kb=0.00

set xrange [0:4]

f(x,a,b)=a/(a+b)*(1-exp(-(a+b)*x))
g1(x)=(kp1/(kp1+kd1))*(1-exp(-(kp1+kd1)*x))
g2(x)=(kp2/(kp2+kd2))*(1-exp(-(kp2+kd2)*x))
g3(x)=(kp3/(kp3+kd3))*(1-exp(-(kp3+kd3)*x))
g4(x)=(kp4/(kp4+kd4))*(1-exp(-(kp4+kd4)*x))
g5(x)=(kp5/(kp5+kd5))*(1-exp(-(kp5+kd5)*x))
g6(x)=(kp6/(kp6+kd6))*(1-exp(-(kp6+kd6)*x))
g7(x)=(kp7/(kp7+kd7))*(1-exp(-(kp7+kd7)*x))

set fit logfile '/dev/null'
set print 'fitPVSATP_PanelE.dat'

fit [0:4] g1(x) './PanelE_101/test1.dat' u 1:2 via kp1,kd1
print 100, kp1, kd1

fit [0:4] g2(x) './PanelE_102/test1.dat' u 1:2 via kp2,kd2
print 86, kp2, kd2

fit [0:4] g3(x) './PanelE_103/test1.dat' u 1:2 via kp3,kd3
print 75, kp3, kd3

fit [0:4] g4(x) './PanelE_104/test1.dat' u 1:2 via kp4,kd4
print 60, kp4, kd4

fit [0:4] g5(x) './PanelE_105/test1.dat' u 1:2 via kp5,kd5
print 50, kp5, kd5

fit [0:4] g6(x) './PanelE_106/test1.dat' u 1:2 via kp6,kd6
print 40, kp6, kd6

fit [0:4] g7(x) './PanelE_107/test1.dat' u 1:2 via kp7,kd7
print 25, kp7, kd7

p g1(x) w l ls 8,\
  g2(x) w l ls 8,\
  g3(x) w l ls 8,\
  g4(x) w l ls 8,\
  g5(x) w l ls 8,\
  g6(x) w l ls 8,\
  g7(x) w l ls 8,\
  './PanelE_101/test1.dat' u 1:COL1 w l ls 1,\
  './PanelE_102/test1.dat' u 1:COL1 w l ls 2,\
  './PanelE_103/test1.dat' u 1:COL1 w l ls 3,\
  './PanelE_104/test1.dat' u 1:COL1 w l ls 4,\
  './PanelE_105/test1.dat' u 1:COL1 w l ls 5,\
  './PanelE_106/test1.dat' u 1:COL1 w l ls 6,\
  './PanelE_107/test1.dat' u 1:COL1 w l ls 7

