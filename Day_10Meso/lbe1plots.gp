set title 'Figure 1: density surface plot'
set view 60, 30, 1, 1
set yrange[0:100]
set xrange[0:100]
set dgrid3d 100,100,1
set contour base
set cntrparam cubicspline
set cntrparam levels 100
set style data lines
unset clabel
set xlabel 'x'
set ylabel 'y'
splot 'data_density.dat' u 1:2:3 nocontour

pause -1 ""


set title 'Figure 2: density contour plot'
set view map
set yrange[0:100]
set xrange[0:100]
set dgrid3d 100,100,1
set contour base
set cntrparam cubicspline
set cntrparam levels 100
set style data lines
unset clabel
set xlabel 'x'
set ylabel 'y'
splot 'data_density.dat' u 1:2:3 w pm3d nocontour, 'data_density.dat' u 1:2:3 with lines lc rgb "black" nosurface notitle

pause -1 ""

set title 'Figure 3: stream function surface plot'
set view 60, 30, 1, 1
set yrange[0:100]
set xrange[0:100]
set dgrid3d 100,100,1
set contour base
set cntrparam cubicspline
set cntrparam levels 100
set style data lines
unset clabel
set xlabel 'x'
set ylabel 'y'
splot 'data_stfn.dat' u 1:2:3 nocontour

pause -1 ""


set title 'Figure 4: streamline (stream function contour plot)'
set view map
set yrange[0:100]
set xrange[0:100]
set dgrid3d 100,100,1
set contour base
set cntrparam cubicspline
set cntrparam levels 100
set style data lines
unset clabel
set xlabel 'x'
set ylabel 'y'
splot 'data_stfn.dat' u 1:2:3 w pm3d nocontour, 'data_stfn.dat' u 1:2:3 with lines lc rgb "black" nosurface notitle

pause -1 ""

