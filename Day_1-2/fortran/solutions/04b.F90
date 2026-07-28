!! Author: Alin M Elena
!! Date: 13-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
module mysum_m
  implicit none
  private
  public :: mySum
  public :: mySum_s

contains

  integer function mySum(n)
    integer, intent(in) :: n

    integer :: i,s
    s = 0
    do i = 1,n
      s = s + i*i
    end do
    mySum = s
  end function mySum

  subroutine mySum_s(s,n)
    integer, intent(in) :: n
    integer, intent(out) :: s

    s = mySum(n)
  end subroutine mySum_s
end module mysum_m

program ex04b
  use mysum_m
  implicit none

  integer :: n,d

  n = 15
  write(*,*) "sum is: ",mySum(n)
  call mySum_s(d,n)
  write(*,*) "sum is: ",d

end program ex04b
