!! Author: Alin M Elena
!! Date: 13-07-2021
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
program ex03c
  implicit none

  integer :: n,i,d

  n = 15
  d = 0
  i = 1
  do while(i<=n)
    if (mod(i,2)==1) then
      i = i + 1
      cycle
     else
      d = d + i*i
      i = i + 1
    end if
  end do
  write(*,*) "sum is: ",d
end program ex03c
