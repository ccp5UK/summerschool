!! Author: Alin M Elena
!! Date: 14-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
module vector_m
  use iso_fortran_env, only : dp=> real64
  implicit none
  private

  type, public :: vector_type
    integer :: n
    real(dp) :: a(100)
  end type 


  public :: norm
  public :: normalize
  public :: positive

contains

  real(dp) function norm(a)

    type(vector_type), intent(in) :: a
    norm = norm2(a%a(1:a%n))

  end function norm 

  subroutine normalize(a)
    type(vector_type), intent(inout) :: a

    a%a(1:a%n) = a%a(1:a%n)/norm(a)

  end subroutine normalize

  subroutine positive(a)
    type(vector_type), intent(inout) :: a

    integer :: i

    do i =1,a%n
      if(a%a(i)<0.0_dp) a%a(i)=-a%a(i)
    end do 
  end subroutine positive

end module vector_m

program test
  use iso_fortran_env, only : dp=> real64
  use vector_m

  implicit none
  integer, parameter :: n=12
  type(vector_type) :: a
  integer :: i

  a%n = n 
  a%a(1:n) = [((i-7.0_dp)/4.0_dp,i=1,n)]
  print *,"a: ",a%a
  print *, "norm: ",norm(a)
  call normalize(a)
  print *,"normalized: ",a%a
  call positive(a)
  print *,"positive: ",a%a
end program test
