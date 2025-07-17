! maths_module.f90
! Some basic self-contained mathematical utility routines for Monte Carlo programs
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module maths_module

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  implicit none
  private

  public :: metropolis, random_integer, random_normal, pbc, trigger, tridiagonal

  ! Generic interface for pbc functions
  interface pbc
     module procedure pbc_one  ! for rank 1 vector (3) of a single position
     module procedure pbc_many ! for rank 2 matrix (3,n) of many positions
  end interface pbc

contains

  function pbc_one ( r, box ) result ( p )
    implicit none
    real, dimension(3), intent(in) :: r   ! Uncorrected position vector
    real, dimension(3), intent(in) :: box ! Simulation box lengths
    real, dimension(3)             :: p   ! Corrected position vector

    ! Applies periodic boundary / minimum image correction to vector r

    p = r - anint ( r / box ) * box

  end function pbc_one

  function pbc_many ( r, box ) result ( p )
    implicit none
    real, dimension(:,:), intent(in)             :: r   ! Uncorrected position vectors (3,n)
    real, dimension(3),   intent(in)             :: box ! Simulation box lengths
    real, dimension(size(r,dim=1),size(r,dim=2)) :: p   ! Corrected position vector (3,n)

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
    real, intent(in) :: delta ! Negative of argument of exponential
    logical          :: acc   ! Returns decision

    ! Set acc=.true. (accept a trial move) with probability min[1,exp(-delta)]
    ! We regard a value delta > 0 as "uphill" while delta < 0 is "downhill"
    
    real            :: zeta
    real, parameter :: exponent_guard = 75.

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

    real :: zeta

    if ( n <= 0 ) stop 'Non-positive argument in random_integer'

    call random_number ( zeta )
    i = 1 + floor ( n*zeta )

    ! Guard against small danger of roundoff
    if ( i < 1 ) i = 1 
    if ( i > n ) i = n

  end function random_integer

  function random_normal ( mean, std ) result ( r )
    implicit none
    real             :: r    ! Returns a normally-distributed random number with
    real, intent(in) :: mean ! specified mean and
    real, intent(in) :: std  ! specified standard deviation

    ! Box-Muller transform produces numbers in pairs
    ! r1 = sqrt(-2.0*log(zeta(1)))*cos(twopi*zeta(2))
    ! r2 = sqrt(-2.0*log(zeta(1)))*sin(twopi*zeta(2))
    ! For this simple application we use one and throw away the other

    real, dimension(2) :: zeta
    real, parameter    :: pi = 4.0*atan(1.0), twopi = 2.0*pi

    call random_number (zeta)                      ! Two uniformly distributed random numbers
    r = sqrt(-2.0*log(zeta(1)))*cos(twopi*zeta(2)) ! Normal, with mean=0, std=1
    r = mean + std * r                             ! Normal, with desired mean, std

  end function random_normal

  function trigger ( x ) result ( yes )
    implicit none
    integer, intent(in) :: x   ! Typically a step counter
    logical             :: yes ! .true. if x equals a trigger value

    real, parameter :: eps = 1.e-10
    real            :: s, ex

    if ( x <= 0 ) stop 'Non-positive argument in trigger'
    
    ! Get x into standard form: s * 10**ex where ex is an integer
    ex = floor(log10(real(x))+tiny(1.0))
    s = x * 10**(-ex)

    ! Trigger values are s = 1, 2, 5 within a small tolerance
    yes = (abs(s-1.)<eps) .or. (abs(s-2.)<eps) .or. (abs(s-5.)<eps)

  end function trigger

  function tridiagonal ( a, b, c, r ) result ( u )
    implicit none
    real, dimension(:), intent(in) :: a ! lower diagonal of matrix M (n-1)
    real, dimension(:), intent(in) :: b ! diagonal of matrix M (n)
    real, dimension(:), intent(in) :: c ! upper diagonal of matrix M (n-1)
    real, dimension(:), intent(in) :: r ! RHS of equations (n)
    real, dimension(size(r))       :: u ! Solution of tridiagonal set (n)

    ! Solves matrix equation M.u = r where M is tridiagonal
    ! Based on the routine tridag_ser in Press, Teukolsky, Vetterling and Flannery
    ! Numerical Recipes in Fortran 2nd edition and Numerical Recipes in Fortran 90
    ! Cambridge University Press (1992 and 1996)
    
    real, dimension(size(r)) :: gam

    integer :: n, j
    real    :: bet

    n = size(r)
    if ( size(a) /= n-1 ) stop 'a dimension error in tridiagonal'
    if ( size(b) /= n   ) stop 'b dimension error in tridiagonal'
    if ( size(c) /= n-1 ) stop 'c dimension error in tridiagonal'

    bet = b(1)
    if (bet == 0.0) stop 'Error 1 in tridiagonal'
    u(1) = r(1) / bet

    ! Decomposition and forward substitution
    do j = 2, n
       gam(j) = c(j-1) / bet
       bet    = b(j) - a(j-1)*gam(j)
       if (bet == 0.0) stop 'Error 2 in tridiagonal'
       u(j) = ( r(j) - a(j-1)*u(j-1) ) / bet
    end do

    ! Back substitution
    do j = n-1, 1, -1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do
    
  end function tridiagonal

end module maths_module
