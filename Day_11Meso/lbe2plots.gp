set title 'Figure 1a: density surface'
set view 60, 30, 1, 1
set size 1,1
set xrange[0:130]
set yrange[0:130]
set dgrid3d 131,131,1
set xlabel 'x'
set ylabel 'y'
unset surface
set pm3d implicit at s
set pm3d scansbackward
splot 'data_density.dat' u 1:2:3 nocontour notitle

pause -1 ""

set title 'Figure 1b: density map'
set view map
set yrange[0:130]
set xrange[0:130]
set dgrid3d 131,131,1
set xlabel 'x'
set ylabel 'y'
set size square
set surface
splot 'data_density.dat' u 1:2:3 w pm3d notitle

pause -1 ""


set title 'Figure 2a: phase index surface'
set view 60, 30, 1, 1
set size 1,1
set xrange[0:130]
set yrange[0:130]
set dgrid3d 131,131,1
set xlabel 'x'
set ylabel 'y'
unset surface
set pm3d implicit at s
set pm3d scansbackward
splot 'data_rhoN.dat' u 1:2:3 nocontour notitle

pause -1 ""


set title 'Figure 2b: phase index map'
set view map
set yrange[0:130]
set xrange[0:130]
set dgrid3d 131,131,1
set xlabel 'x'
set ylabel 'y'
set surface
splot 'data_rhoN.dat' u 1:2:3 w pm3d nocontour notitle

pause -1 ""


set title 'Figure 3a: curvature surface'
set view 60, 30, 1, 1
set size 1,1
set xrange[0:130]
set yrange[0:130]
set dgrid3d 131,131,1
set xlabel 'x'
set ylabel 'y'
unset surface
set pm3d implicit at s
set pm3d scansbackward
splot 'data_nK.dat' u 1:2:3 nocontour notitle

pause -1 ""


set title 'Figure 3b: curvature map'
set view map
set yrange[0:130]
set xrange[0:130]
set dgrid3d 131,131,1
set xlabel 'x'
set ylabel 'y'
set surface
splot 'data_nK.dat' u 1:2:3 w pm3d

pause -1 ""

set title 'Figure 4a: velocity modulus surface'
set view 60, 30, 1, 1
set size 1,1
set xrange[0:130]
set yrange[0:130]
set dgrid3d 131,131,1
set xlabel 'x'
set ylabel 'y'
unset surface
set pm3d implicit at s
set pm3d scansbackward
splot 'data_psi.dat' u 1:2:3 nocontour notitle

pause -1 ""


set title 'Figure 4b: velocity modulus map'
set view map
set yrange[0:130]
set xrange[0:130]
set dgrid3d 131,131,1
set xlabel 'x'
set ylabel 'y'
set surface
splot 'data_psi.dat' u 1:2:3 w pm3d

pause -1 ""

