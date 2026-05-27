!! Author: Alin M Elena
!! Date: 10-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
program ex02
  use iso_fortran_env, only : dp=> real64
  implicit none

  real(dp) :: a,b,c
  real(dp) :: d

  !a = 1.0_dp;b=-4.0_dp;c=3.0_dp
  !a = 0.0_dp;b=-4.0_dp;c=3.0_dp
  !a = 1.0_dp;b=-2.0_dp;c=1.0_dp
  a = 1.0_dp;b=2.0_dp;c=2.0_dp

  d = b*b - 4.0_dp*a*c
  print *,"discriminant ",d
  if (abs(a)<epsilon(1.0_dp)) then
    write(*,*) "One root, linear equation: ", -c/b
  elseif(abs(d)<epsilon(1.0_dp)) then
    write(*,*) "One root: ", -b/(2.0_dp*a)
  elseif(d<tiny(1.0_dp)) then
    write(*,*) "two complex roots: ", (-b+sqrt(cmplx(0.0_dp,d)))/(2.0_dp*a),&
      (-b-sqrt(cmplx(0.0_dp,d)))/(2.0_dp*a)
  else
    write(*,*) "two real roots: ", (-b+sqrt(d))/(2.0_dp*a),&
      (-b-sqrt(d))/(2.0_dp*a)

  end if

end program ex02
