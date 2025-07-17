!! Author: Alin M Elena
!! Date: 13-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
program ex03c
  implicit none

  integer :: n,i,d

  n = 15
  d = 0
  i = 2
  do while(i<=n)
    d = d + i*i
    i = i + 2
  end do
  write(*,*) "sum is: ",d
end program ex03c
