set term win
set out
set style line 1 lw 3
set style line 2 lw 3
set style line 3 lw 3
set style line 4 lw 3
set style line 5 lw 3
set style line 6 lw 3
set style line 7 lw 3
set ylab ''
set xlab ''
set nokey
set size square
#set noxtics
#set noytics
#set notics
set key l w 2
set pm3d map
unset surf
set contour
set cntrparam levels discrete .5
set terminal png enhanced font Arial 26 size 1028,1028 crop
set out '3dprf_10keVc_Ex2.png'
splot 'mantis_demo_CsI_HR_10keV3det.out' t '10 keV' w l lt 0 lw 2

set surf
unset contour
set pm3d map
set out '3dprf_10keV_Ex2.png'
splot 'mantis_demo_CsI_HR_10keV3det.out' t '10 keV' #w l lt 0 lw 2

set out
