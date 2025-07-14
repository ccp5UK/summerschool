! mc_muca_module.f90
! Routine to carry out insertion/deletion moves for muca program
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module mc_muca_module

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! ACKNOWLEDGEMENT
  ! An earlier code written by Leo Lue, on which this was based, is gratefully acknowledged.
  ! Needless to say, he is in no way responsible for any mistakes in this code.

  implicit none
  private

  public :: n_move

contains

  subroutine n_move ( beta, phi, box, n, r, u, w, accepted ) ! Particle insertion or deletion move
    use maths_module,     only : random_integer, metropolis
    use potential_module, only : pot_one
    implicit none

    real,                    intent(in)    :: beta     ! Inverse temperature
    real,    dimension(0:),  intent(in)    :: phi      ! Array of weights (0:nmax)
    real,    dimension(3),   intent(in)    :: box      ! Box lengths
    integer,                 intent(inout) :: n        ! Current number of atoms
    real,    dimension(:,:), intent(inout) :: r        ! Array of coordinates (3,nmax)
    real,                    intent(inout) :: u        ! Total potential energy
    real,                    intent(inout) :: w        ! Total virial
    integer,                 intent(inout) :: accepted ! Accepted move counter

    ! Carries out single-atom insertion or deletion
    ! The ensemble is defined by beta and weights phi(n)
    ! The maximum number of atoms, nmax, is defined by the dimensions of phi and r
    ! The actual number of atoms, n, must not exceed nmax
    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values.
    ! If the move is rejected, they will be returned unchanged.

    real, dimension(3) :: r_new ! Temporary position array

    integer :: nmax, i
    real    :: zeta, v
    real    :: u_new, w_new, u_old, w_old, delta
    logical :: overlap

    if (   size ( r, dim=1 ) /= 3 ) stop 'r dimension error in n_move'
    nmax = size ( r, dim=2 )
    if ( ubound(phi,dim=1) /= nmax ) stop 'phi dimension error in n_move'
    if ( n > nmax  ) stop 'n error in n_move'
    v  = product ( box ) ! Volume

    call random_number ( zeta ) ! Choose between insertion and deletion

    if ( zeta > 0.5 ) then ! Attempt a particle insertion move

       if ( n == nmax ) return ! Reject: can't exceed nmax

       ! Choose random position
       call random_number ( r_new )
       r_new = box * ( r_new - 0.5 )

       ! Check for overlap, and calculate energy and virial, for interaction with other atoms
       call pot_one ( box, r(:,:n), r_new, 0, u_new, w_new, overlap )
       if ( overlap ) return ! Reject move on overlap

       delta = beta  *u_new                         ! Delta U contribution
       delta = delta - log (v/real(n+1))            ! Contribution related to (N+1)/V
       delta = delta - beta * ( phi(n+1) - phi(n) ) ! Delta Phi contribution

       if ( metropolis ( delta ) ) then ! Accept Metropolis test
          n        = n + 1        ! Increase number of atoms
          r(:,n)   = r_new        ! Add new atom at the end
          u        = u + u_new    ! Add new interaction potential
          w        = w + w_new    ! Add new interaction virial
          accepted = accepted + 1 ! Increment move counter
       end if

    else ! Attempt a particle deletion move

       if ( n == 0 ) return ! Reject: can't go below 0

       ! Choose target atom at random
       i = random_integer ( n )

       ! Compute current interactions
       call pot_one ( box, r(:,:n), r(:,i), i, u_old, w_old, overlap )
       if ( overlap ) stop 'Overlap in current configuration' ! Should never happen

       delta = -beta * u_old                        ! Delta U contribution
       delta = delta + log (v/real(n))              ! Contribution related to N/V
       delta = delta - beta * ( phi(n-1) - phi(n) ) ! Delta Phi contribution

       if ( metropolis ( delta ) ) then ! Accept Metropolis test
          r(:,i)   = r(:,n)       ! Overwrite atom i with atom n
          n        = n - 1        ! Decrease number of atoms
          u        = u - u_old    ! Subtract old interaction potential
          w        = w - w_old    ! Subtract old interaction virial
          accepted = accepted + 1 ! Increment move counter
       end if

    end if

  end subroutine n_move

end module mc_muca_module
