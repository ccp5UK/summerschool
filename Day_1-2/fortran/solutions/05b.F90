!! Author: Alin M Elena
!! Date: 13-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
program ex05
  use iso_fortran_env, only : dp=> real64
  implicit none
  integer :: n
  real(dp) :: eps,diff,old,new
  n = 2
  eps = 1.0e-6
  new=simp(n,0.0_dp,2.0_dp)
  diff =huge(1.0_dp)
  do while (diff>eps)
    n = n + 2
    old = new
    new = simp(n,0.0_dp,2.0_dp)
    diff = abs(new-old)
  end do
  print *, "integral value ", new, "with ",n

contains
  real(dp) function f(x)
    real(dp),intent(in) :: x

    f = x*x -2.0_dp*x +1.0_dp
!f = x*x*x -2.0_dp*exp(-x*x)+1.0_dp
  end function f

  real(dp) function simp(n,a,b)
    integer, intent(in) :: n
    real(dp),intent(in) :: a,b

    real(dp) :: h,s
    integer :: i

    h = (b-a)/real(n,dp)

    s =0.0_dp
    do i = 1, n/2
      s = s + f(a+h*(2*i-2))+4.0_dp*f(a+h*(2*i-1))+f(a+h*2*i)
    end do
    simp =h/3.0_dp * s

  end function simp
end program ex05
