!! Author: Alin M Elena
!! Date: 13-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
program ex03a
  implicit none

  integer :: n,i,d

  n = 15
  d = 0
  do i = 1,n
    if ( mod(i,2) == 1 ) cycle
    d = d + i*i
  end do
  write(*,*) "sum is: ",d
end program ex03a
