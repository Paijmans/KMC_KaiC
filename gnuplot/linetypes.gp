# Linetype library:
# 1 - solid line __________
# 2 - dashed as __ __ __ __
# 3 - dashed as _ _ _ _ _ _
# 4 - dashed as ...........
# 5 - dashed as __.__.__.__
# 6 - dashed as _._._._._._


set style line 11 lw 5 lt 1 lc rgb '#FF0000'
set style line 12 lw 5 lt 1 lc rgb '#00FF00'
set style line 13 lw 5 lt 1 lc rgb '#0000FF'
set style line 14 lw 5 lt 1 lc rgb '#FFD700'
set style line 15 lw 5 lt 1 lc rgb '#000000'


# Phosphorylation level to ATP/ADP level traces by Rust
set style line 21 lw 5 lt 1 lc rgb '#ed1a93'
set style line 22 lw 5 lt 1 lc rgb '#f58612'
set style line 23 lw 5 lt 1 lc rgb '#879124'
set style line 24 lw 5 lt 1 lc rgb '#26913f'
set style line 25 lw 5 lt 1 lc rgb '#3f50ab'
set style line 26 lw 5 lt 1 lc rgb '#7c2b98'
set style line 27 lw 5 lt 1 lc rgb '#231f20'
set style line 28 lw 5 lt 1 lc rgb '#bebe2c'
set style line 29 lw 5 lt 1 lc rgb '#6cbe2c'




#Gray lines
set style line 111 lw 4 lt 1 lc rgb '#000000'
set style line 112 lw 4 lt 2 lc rgb '#000000'
set style line 113 lw 4 lt 3 lc rgb '#000000'
set style line 114 lw 4 lt 4 lc rgb '#000000'
set style line 115 lw 5 lt 5 lc rgb '#000000'

set style line 211 lw 1 lt 2 lc rgb '#000000'
set style line 212 lw 2 lt 2 lc rgb 'gray90'
set style line 213 lw 2 lt 2 lc rgb 'gray70'

set style line 214 lw 2 lt 1 lc rgb 'gray30'
set style line 224 lw 2 lt 2 lc rgb 'gray30'


#Arrow styles
set style arrow 1  head filled size screen 0.025,30,45 lw 1 lt 1 lc rgbcolor 'black'
set style arrow 2  head nofilled size screen 0.03,15 lw 1 lt 1 lc rgbcolor 'black'
set style arrow 3  head filled size screen 0.03,15,45 lw 1 lt 1 lc rgbcolor 'black'
set style arrow 31 head filled size screen 0.03,15,45 lw 1 lt 1 lc rgbcolor 'white'
set style arrow 32  head filled size screen 0.02,45,45 lw 1 lt 1 lc rgbcolor 'black'
set style arrow 4  heads size screen 0.008,90 lw 2 lt 1 lc rgbcolor 'black'
set style arrow 5  heads filled size screen 0.03,15,45 lw 2 lt 1 lc rgbcolor 'green'
set style arrow 6  heads filled size screen 0.03,15,45 lw 2 lt 2 lc rgbcolor 'red'
set style arrow 7  nohead ls 214
set style arrow 8  nohead ls 224
set style arrow 9  head filled size screen 0.03,15,45 ls 214
set style arrow 10 head filled size screen 0.03,15,45 ls 224

