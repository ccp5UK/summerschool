! maths_module.f90
! Some basic self-contained mathematical utility routines
module maths_module

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env, only : dp => real64

  implicit none

  ! Public routines
  private
  public :: metropolis, random_integer, random_normal, pbc, trigger, tridiagonal, cbrt, cumsum, cross_products

  ! Generic interface for pbc functions
  interface pbc
     module procedure pbc_one  ! for rank 1 vector (3) of a single position
     module procedure pbc_many ! for rank 2 matrix (3,n) of many positions
  end interface pbc

contains

  function pbc_one ( r, box ) result ( p )
    implicit none
    real(dp), intent(in) :: r(3)   ! Uncorrected position vector
    real(dp), intent(in) :: box(3) ! Simulation box lengths

    real(dp) :: p(3) ! Corrected position vector

    ! Applies periodic boundary / minimum image correction to vector r

    p = r - anint ( r / box ) * box

  end function pbc_one

  function pbc_many ( r, box ) result ( p )
    implicit none
    real(dp), intent(in) :: r(:,:) ! Uncorrected position vectors (3,n)
    real(dp), intent(in) :: box(3) ! Simulation box lengths

    real(dp) :: p(size(r,dim=1),size(r,dim=2)) ! Corrected position vector (3,n)

    ! Applies periodic boundary / minimum image correction to n vectors in r(3,n)

    integer :: n, i

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in pbc_many'
    n = size(r,dim=2)

    ! We could use this statement instead of the following loop
    ! p = r - anint ( r / spread(box,dim=2,ncopies=n) ) * spread(box,dim=2,ncopies=n)
    do i = 1, n
       p(:,i) = r(:,i) - anint ( r(:,i) / box ) * box
    end do

  end function pbc_many
  
  function metropolis ( delta ) result ( acc ) ! Conduct Metropolis test, with safeguards
    implicit none
    real(dp), intent(in) :: delta ! Negative of argument of exponential
    logical              :: acc   ! Returns decision

    ! Set acc=.true. (accept a trial move) with probability min[1,exp(-delta)]
    ! We regard a value delta > 0 as "uphill" while delta < 0 is "downhill"
    
    real(dp)            :: zeta
    real(dp), parameter :: exponent_guard = 75.0_dp

    if ( delta > exponent_guard ) then ! Too high, reject without evaluating
       acc = .false.
    else if ( delta < 0 ) then ! Downhill, accept without evaluating
       acc = .true.
    else
       call random_number ( zeta ) ! Uniform random number in range (0,1)
       acc = exp(-delta) > zeta ! Metropolis test
    end if

  end function metropolis

  function random_integer ( n ) result ( i )
    implicit none
    integer             :: i ! Returns a uniformly distributed random integer
    integer, intent(in) :: n ! in the range [1,n] inclusive

    real(dp) :: zeta

    if ( n <= 0 ) stop 'Non-positive argument in random_integer'

    call random_number ( zeta )
    i = 1 + floor ( n*zeta )

    ! Guard against small danger of roundoff
    if ( i < 1 ) i = 1 
    if ( i > n ) i = n

  end function random_integer

  function random_normal ( mean, std ) result ( r )
    implicit none
    real(dp)             :: r    ! Returns a normally-distributed random number with
    real(dp), intent(in) :: mean ! specified mean and
    real(dp), intent(in) :: std  ! specified standard deviation

    ! Box-Muller transform produces numbers in pairs
    ! r1 = sqrt(-2.0*log(zeta(1)))*cos(twopi*zeta(2))
    ! r2 = sqrt(-2.0*log(zeta(1)))*sin(twopi*zeta(2))
    ! For this simple application we use one and throw away the other

    real(dp) :: zeta(2)

    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp), twopi = 2.0_dp*pi

    call random_number (zeta)                    ! Two uniformly distributed random numbers
    r = sqrt(-2*log(zeta(1)))*cos(twopi*zeta(2)) ! Normal, with mean=0, std=1
    r = mean + std * r                           ! Normal, with desired mean, std

  end function random_normal

  function trigger ( x ) result ( yes )
    implicit none
    integer, intent(in) :: x   ! Typically a step counter
    logical             :: yes ! .true. if x equals a trigger value

    real(dp), parameter :: eps = 1.0e-10_dp
    real(dp)            :: s, ex

    if ( x <= 0 ) stop 'Non-positive argument in trigger'
    
    ! Get x into standard form: s * 10**ex where ex is an integer
    ex = floor(log10(real(x))+tiny(1.0_dp))
    s  = x * 10**(-ex)

    ! Trigger values are s = 1, 2, 5 within a small tolerance
    yes = (abs(s-1.0_dp)<eps) .or. (abs(s-2.0_dp)<eps) .or. (abs(s-5.0_dp)<eps)

  end function trigger

  function tridiagonal ( a, b, c, r ) result ( u )
    implicit none
    real(dp), intent(in) :: a(:) ! lower diagonal of matrix M (n-1)
    real(dp), intent(in) :: b(:) ! diagonal of matrix M (n)
    real(dp), intent(in) :: c(:) ! upper diagonal of matrix M (n-1)
    real(dp), intent(in) :: r(:) ! RHS of equations (n)

    real(dp) :: u(size(r)) ! Solution of tridiagonal set (n)

    ! Solves matrix equation M.u = r where M is tridiagonal
    ! Based on the routine tridag_ser in Press, Teukolsky, Vetterling and Flannery
    ! Numerical Recipes in Fortran 2nd edition and Numerical Recipes in Fortran 90
    ! Cambridge University Press (1992 and 1996)
    
    real(dp) :: gam(size(r))

    integer  :: n, j
    real(dp) :: bet

    n = size(r)
    if ( size(a) /= n-1 ) stop 'a dimension error in tridiagonal'
    if ( size(b) /= n   ) stop 'b dimension error in tridiagonal'
    if ( size(c) /= n-1 ) stop 'c dimension error in tridiagonal'

    bet = b(1)
    if (bet == 0.0_dp) stop 'Error 1 in tridiagonal'
    u(1) = r(1) / bet

    ! Decomposition and forward substitution
    do j = 2, n
       gam(j) = c(j-1) / bet
       bet    = b(j) - a(j-1)*gam(j)
       if (bet == 0.0_dp) stop 'Error 2 in tridiagonal'
       u(j) = ( r(j) - a(j-1)*u(j-1) ) / bet
    end do

    ! Back substitution
    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do
    
  end function tridiagonal
  
  elemental function cbrt ( x ) result ( y ) ! Compute cube root of argument
    implicit none
    real(dp), intent(in) :: x ! Argument
    real(dp)             :: y ! Function return y = x**(1/3)

    real(dp), parameter :: one_third = 1.0_dp / 3.0_dp

    ! Computes cube root of real argument, allowing for negative values
    ! Alternatively could use sign(exp(log(abs(x))/3),x)
    
    y = sign(abs(x)**one_third,x)

  end function cbrt

  function cumsum ( x ) result ( c ) ! Cumulative sum function for supplied array
    real(dp), intent(in) :: x(:)

    real(dp) :: c(size(x))
    integer  :: j

    if ( size(x) == 0 ) return

    c(1) = x(1)
    do j = 2, size(x)
       c(j) = c(j-1) + x(j)
    end do

  end function cumsum

  function cross_products ( a, b ) result ( c )
    implicit none
    real(dp), intent(in) :: a(:,:), b(:,:) ! Two supplied sets of vectors (3,n)

    real(dp) :: c(size(a,dim=1),size(a,dim=2)) ! Returns vector cross products (3,n)

    integer :: i, n

    if ( size(a,dim=1) /= 3 ) stop 'a dimension error in cross_product_many'
    if ( size(b,dim=1) /= 3 ) stop 'b dimension error in cross_product_many'
    n = size(a,dim=2)
    if ( size(b,dim=2) /= n ) stop 'ab dimension error in cross_product_many'

    do i = 1, n
       c(1,i) = a(2,i)*b(3,i) - a(3,i)*b(2,i)
       c(2,i) = a(3,i)*b(1,i) - a(1,i)*b(3,i)
       c(3,i) = a(1,i)*b(2,i) - a(2,i)*b(1,i)
    end do
  end function cross_products

end module maths_module
