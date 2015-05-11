reset 
clear
set term pos  enh eps color 24
set key left
set auto
set ylab 'Frequency'
set xlab 'Number of optical quanta per primary x ray'

set out 'phs.eps'
plot   '../demo/mantis_demo_CsI_HRcollect.out' t 'HR' w l lt 1 lw 3  \
, '../demo/mantis_demo_CsI_HLcollect.out' t 'HL' w l lt 2 lw 3 

set auto
set pm3d map corners2color c1
set pal gray gamma 2
unset log
set ylab ''
set xlab ''
set vi 0,0
unset colorbox
unset title
unset key
set size square
set size ratio -1
unset xtics
unset ytics
unset tics
unset surf
set contour
set cntrparam levels discrete .1
set log z
unset log
set term png size 1000,1000 crop

set out '3dprf_HR.png'
splot '../demo/mantis_demo_CsI_HR3det.out' w l  

set out '3dprf_HL.png'
splot '../demo/mantis_demo_CsI_HL3det.out' w l  

unset pm3d
set out
