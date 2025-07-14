! mc_module.f90
! Routines to carry out atom translation moves and test-particle insertion for Monte Carlo programs
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module mc_module

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  implicit none
  private

  public :: r_move, insert

contains

  subroutine r_move ( beta, dr_max, box, r, u, w, accepted ) ! Single-atom displacement
    use maths_module,     only : random_integer, pbc, metropolis
    use potential_module, only : pot_one
    implicit none

    real,                    intent(in)    :: beta     ! Inverse temperature
    real,                    intent(in)    :: dr_max   ! Maximum displacement
    real,    dimension(3),   intent(in)    :: box      ! Box lengths
    real,    dimension(:,:), intent(inout) :: r        ! Array of coordinates (3,n)
    real,                    intent(inout) :: u        ! Total potential energy
    real,                    intent(inout) :: w        ! Total virial
    integer,                 intent(inout) :: accepted ! Accepted move counter

    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values.
    ! If the move is rejected, they will be returned unchanged.

    real, dimension(3) :: r_new, dr

    integer :: i, n
    real    :: u_new, w_new, u_old, w_old, delta
    logical :: overlap

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in r_move'
    n  = size(r,dim=2) ! Number of atoms

    if ( n <= 0 ) return ! Nothing to move

    i = random_integer ( n ) ! Choose atom randomly

    call pot_one ( box, r, r(:,i), i, u_old, w_old, overlap )
    if ( overlap ) stop 'Overlap in current configuration' ! Should never happen

    call random_number ( dr )   ! Choose random displacement
    dr  = ( 2*dr - 1 ) * dr_max ! in range (-dr_max,+dr_max)
    r_new = r(:,i) + dr         ! Trial move to new position 
    r_new = pbc ( r_new, box )  ! Periodic boundary correction

    call pot_one ( box, r, r_new, i, u_new, w_new, overlap )
    if ( overlap ) return ! Reject trial move on overlap

    delta = beta * ( u_new - u_old )

    if ( metropolis ( delta ) ) then ! Accept Metropolis test
       r(:,i)   = r_new                ! Update position
       u        = u + u_new - u_old    ! Update potential
       w        = w + w_new - w_old    ! Update virial
       accepted = accepted + 1         ! Increment move counter
    end if ! End accept Metropolis test

  end subroutine r_move

  function insert ( beta, box, r ) result ( zin ) ! Test particle insertion
    use potential_module, only : pot_one
    implicit none

    real,                    intent(in) :: beta ! Inverse temperature
    real,    dimension(3),   intent(in) :: box  ! Box lengths
    real,    dimension(:,:), intent(in) :: r    ! Array of coordinates (3,n)
    real                                :: zin  ! Widom estimate of 1/z

    integer              :: n
    real,   dimension(3) :: r_new
    real                 :: v, u_new, w_new
    logical              :: overlap

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in insert'
    n  = size(r,dim=2) ! Number of atoms

    zin = 0

    call random_number ( r_new )                             ! 3 random numbers in range (0,1)
    r_new = ( r_new - 0.5 ) * box                            ! Random position in box
    call pot_one ( box, r, r_new, 0, u_new, w_new, overlap ) ! Interaction with real atoms
    if ( overlap ) return ! zin=0 on overlap

    v   = product ( box ) ! Volume
    zin = exp(-beta*u_new) * v / (n+1)

  end function insert

end module mc_module
