!! Author: Alin M Elena
!! Date: 14-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
program pi
  use iso_fortran_env, only : dp=> real64
  implicit none

  integer(kind=4), allocatable :: seed(:)
  real(kind=8) :: r(2)
  integer :: is,i,p
  integer :: n ,c

  is = 13
  call random_seed(size=p)
  allocate(seed(p))
  seed = 17*[(i-is,i=1,p)]
  call random_seed(put=seed)
  deallocate(seed)

  n = 20000000
  c = 0

  do i=1,n
    call random_number(r)
    if (r(1)**2+r(2)**2<1.0_dp) c = c + 1
  end do 
  print *, 4.0_dp*real(c,dp)/real(n,dp)
end program pi
