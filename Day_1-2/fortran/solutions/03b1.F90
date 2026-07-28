!! Author: Alin M Elena
!! Date: 13-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
program ex03b
  implicit none

  integer :: n,i,d

  n = 15
  d = 0
  i = 2
  do
    d = d + i*i
    i = i + 2
    if (i>n) exit
  end do
  write(*,*) "sum is: ",d
end program ex03b
