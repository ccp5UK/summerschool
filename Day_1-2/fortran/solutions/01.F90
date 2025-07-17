!! Author: Alin M Elena
!! Date: 10-07-2019
!! License: GPL-3.0 https://opensource.org/licenses/GPL-3.0
program ex01
  use iso_fortran_env, only : dp=> real64, sp=>real32
  implicit none


  real(kind=dp) :: a
  real(kind=sp) :: b
  complex(kind=dp) :: c 
  character(len=42) :: k
  integer :: h


  a=30.0_dp
  b=5.0_sp
  h = 42
  c=cmplx(a,a/2.0_dp)
  k = "   john smith  "
  write(*,*)"len: ", len(k)
  write(*,*)"trim: ", trim(k), len(trim(k))
  write(*,*)"trim: ", adjustl(k), len(adjustl(k))
  write(*,*) "sin: ",sin(a)," cos: ",cos(a)
  write(*,*)"complex ", c
  write(*,*)"real part", real(c)
  write(*,*)"imaginary part", aimag(c)
  write(*,*)"mod ",mod(h,5),h/5 
  write(*,*)"single precision:"
  write(*,*)"tiny: ",tiny(b)
  write(*,*)"huge: ",huge(b)
  write(*,*)"espsilon: ",epsilon(b)
  write(*,*)"double precision:"
  write(*,*)"tiny: ",tiny(a)
  write(*,*)"huge: ",huge(a)
  write(*,*)"espsilon: ",epsilon(a)
end program 
