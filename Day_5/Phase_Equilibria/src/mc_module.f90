! mc_module.f90
! Routines to carry out moves, swaps, and test-particle insertion for Monte Carlo programs
module mc_module

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! In the "Monte Carlo at Constant Pressure" workshop, you were asked to fix
  ! the v_move routine which carries out volume moves.
  ! The calculation of delta has been FIXED in this file.

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env, only : dp => real64

  implicit none

  ! Public routines
  private
  public :: r_move, v_move, n_swap, v_swap, n_move, insert

contains

  subroutine r_move ( beta, dr_max, box, r, u, w, accepted ) ! Single-atom displacement
    use maths_module,     only : random_integer, pbc, metropolis
    use potential_module, only : pot_one
    implicit none

    real(dp), intent(in)    :: beta     ! Inverse temperature
    real(dp), intent(in)    :: dr_max   ! Maximum displacement
    real(dp), intent(in)    :: box(3)   ! Box lengths
    real(dp), intent(inout) :: r(:,:)   ! Array of coordinates (3,n)
    real(dp), intent(inout) :: u        ! Total potential energy
    real(dp), intent(inout) :: w        ! Total virial
    integer,  intent(inout) :: accepted ! Accepted move counter

    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values.
    ! If the move is rejected, they will be returned unchanged.

    real(dp) :: r_new(3), dr(3)

    integer  :: i, n
    real(dp) :: u_new, w_new, u_old, w_old, delta
    logical  :: overlap

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in r_move'
    n  = size(r,dim=2) ! Number of atoms

    if ( n <= 0 ) return ! Nothing to move

    i = random_integer ( n ) ! Choose atom randomly

    call pot_one ( box, r, r(:,i), i, u_old, w_old, overlap )
    if ( overlap ) stop 'Overlap in current configuration' ! Should never happen

    call random_number ( dr )  ! Random vector with components in range (0,1)
    dr  = (2*dr-1) * dr_max    ! Components now in range (-dr_max,+dr_max)
    r_new = r(:,i) + dr        ! Trial move to new position 
    r_new = pbc ( r_new, box ) ! Periodic boundary correction

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

    real(dp), intent(in) :: beta   ! Inverse temperature
    real(dp), intent(in) :: box(3) ! Box lengths
    real(dp), intent(in) :: r(:,:) ! Array of coordinates (3,n)
    real(dp)             :: zin    ! Widom estimate of 1/z

    integer  :: n
    real(dp) :: r_new(3)
    real(dp) :: v, u_new, w_new
    logical  :: overlap

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in insert'
    n  = size(r,dim=2) ! Number of atoms

    zin = 0

    call random_number ( r_new )                             ! Random vector, components in range (0,1)
    r_new = (2*r_new-1) * box/2                              ! Components now in range (-box(:)/2,box(:)/2)
    call pot_one ( box, r, r_new, 0, u_new, w_new, overlap ) ! Interaction with real atoms
    if ( overlap ) return ! zin=0 on overlap

    v   = product ( box ) ! Volume
    zin = exp(-beta*u_new) * v / (n+1)

  end function insert

  subroutine v_move ( beta, pressure, dv_max, box, r, u, w, accepted ) ! Volume move
    use maths_module,     only : metropolis, cbrt
    use potential_module, only : pot_all
    implicit none

    real(dp), intent(in)    :: beta     ! Inverse temperature
    real(dp), intent(in)    :: pressure ! Pressure
    real(dp), intent(in)    :: dv_max   ! Maximum volume change
    real(dp), intent(inout) :: box(3)   ! Box lengths
    real(dp), intent(inout) :: r(:,:)   ! Array of coordinates (3,n)
    real(dp), intent(inout) :: u        ! Total potential energy
    real(dp), intent(inout) :: w        ! Total virial
    integer,  intent(inout) :: accepted ! Accepted move counter

    ! Volume move in constant-NPT Monte Carlo

    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values for the scaled system.
    ! If the move is rejected, they will be returned unchanged.
    
    real(dp) :: r_new(size(r,dim=1),size(r,dim=2)) ! New (scaled) coordinates

    integer  :: n
    real(dp) :: dv, delta
    real(dp) :: box_new(3)
    real(dp) :: v_new, v, r_scale, v_scale
    real(dp) :: u_new, w_new
    logical  :: overlap

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in v_move'
    n  = size(r,dim=2) ! Number of atoms

    call random_number ( dv )    ! Randoom number uniform in range (0,1)
    dv      = (2*dv-1) * dv_max  ! Now uniform in range (-dv_max,+dv_max)
    v       = product ( box )    ! Current (or old) volume
    v_new   = v + dv             ! New volume
    v_scale = v_new / v          ! Scaling factor for volume
    r_scale = cbrt(v_scale)      ! Scaling factor for coordinates and box lengths
    box_new = box * r_scale      ! New box lengths
    r_new   = r   * r_scale      ! New coordinates

    call pot_all ( box_new, r_new, u_new, w_new, overlap ) 
    if ( overlap ) return ! Reject move on overlap

    delta = beta * ( u_new - u )         ! This is the exp(-beta*Delta U) part
    delta = delta + beta * pressure * dv ! This is the exp(-beta*P*Delta V) part (corrected)
    delta = delta - n * log ( v_scale )  ! This is the (Vnew/Vold)**N part (no correction needed)

    if ( metropolis ( delta ) ) then ! Accept Metropolis test
       u        = u_new        ! Update potential
       w        = w_new        ! Update virial
       box      = box_new      ! Update box lengths
       r        = r_new        ! Update coordinates
       accepted = accepted + 1 ! Increment move counter
    end if

  end subroutine v_move

  subroutine n_swap ( beta, box, n, r, u, w, accepted ) ! Particle exchange move
    use maths_module,     only : random_integer, metropolis
    use potential_module, only : pot_one
    implicit none

    real(dp), intent(in)    :: beta     ! Inverse temperature
    real(dp), intent(in)    :: box(3,2) ! Box lengths in both systems
    integer,  intent(inout) :: n(2)     ! Number of atoms in both systems
    real(dp), intent(inout) :: r(:,:)   ! Combined array of coordinates (3,n(1)+n(2))
    real(dp), intent(inout) :: u(2)     ! Total potential energy in both systems
    real(dp), intent(inout) :: w(2)     ! Total virial in both systems
    integer,  intent(inout) :: accepted ! Accepted move counter

    ! Particle exchange move in Gibbs ensemble Monte Carlo
    
    ! Attempts to carry out a particle exchange move, with the direction 1->2 or 2->1 picked randomly.
    ! Both systems are contained within the r array: r(1:n(1)) and r(n(1)+1,n(1)+n(2)).
    ! The particle being moved is selected randomly within the origin box.
    ! The insertion position is chosen at random in the destination box.
    ! Note that we could take advantage of this to do a Widom test particle calculation at the same time.
    ! However, in this program, we are doing test particle insertions in a separate routine,
    ! which is less efficient, but more clear.

    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values.
    ! If the move is rejected, they will be returned unchanged.

    ! NB variables such as box, n, u, w have an extra dimension, corresponding to the two simulated systems.

    integer  :: i, k, ii, choice
    real(dp) :: r_new(3)
    real(dp) :: v(2)
    real(dp) :: u_new, w_new, u_old, w_old, delta
    logical  :: overlap

    if ( any ( shape(r) /= [3,sum(n)] ) ) stop 'r dimension error in n_swap'

    v = product ( box, dim=1 ) ! Both box volumes

    choice = random_integer(2) ! Choose direction of swap

    if ( choice == 1 ) then ! Try swapping 1->2

       if ( n(1) <= 1 ) return ! Disallow n(1)->0

       i = random_integer ( n(1) ) ! Choose atom in system 1 at random
       call pot_one ( box(:,1), r(:,:n(1)), r(:,i), i, u_old, w_old, overlap )
       if ( overlap ) stop 'Overlap in current configuration' ! Should never happen
       call random_number ( r_new )     ! Random vector with components in range (0,1)
       r_new = (2*r_new-1) * box(:,2)/2 ! Components now in range of system 2 box
       call pot_one ( box(:,2), r(:,n(1)+1:), r_new, 0, u_new, w_new, overlap )
       if ( overlap ) return ! Reject move on overlap

       delta = beta * ( u_new - u_old )        ! Delta U contribution
       delta = delta - log ( v(2) / (n(2)+1) ) ! Contribution from creation in 2
       delta = delta + log ( v(1) / n(1)     ) ! Contribution from destruction in 1

       if ( metropolis ( delta ) ) then ! Accept Metropolis test
          k        = n(1)         ! k is last atom in system 1
          r(:,i)   = r(:,k)       ! Replace atom i by atom k
          r(:,k)   = r_new        ! Put new particle in position k
          n        = n + [-1,1]   ! Move boundary down so k is now in system 2
          u(1)     = u(1) - u_old ! Update potential in system 1
          w(1)     = w(1) - w_old ! Update virial in system 1
          u(2)     = u(2) + u_new ! Update potential in system 2
          w(2)     = w(2) + w_new ! Update virial in system 2
          accepted = accepted + 1 ! Increment move counter
       end if

    else ! Try swapping 2->1

       if ( n(2) <= 1 ) return ! Disallow n(2)->0

       i  = random_integer ( n(2) ) ! Choose atom in system 2 at random
       ii = i + n(1)                ! This is its position in the r array
       call pot_one ( box(:,2), r(:,n(1)+1:), r(:,ii), i, u_old, w_old, overlap )
       if ( overlap ) stop 'Overlap in current configuration' ! Should never happen
       call random_number ( r_new )     ! Random vector with components in range (0,1)
       r_new = (2*r_new-1) * box(:,1)/2 ! Components now in range of system 1 box
       call pot_one ( box(:,1), r(:,:n(1)), r_new, 0, u_new, w_new, overlap )
       if ( overlap ) return ! Reject move on overlap

       delta = beta * ( u_new - u_old )        ! Delta U contribution
       delta = delta - log ( v(1) / (n(1)+1) ) ! Contribution from creation in 1
       delta = delta + log ( v(2) / n(2)     ) ! Contribution from destruction in 2

       if ( metropolis ( delta ) ) then ! Accept Metropolis test
          k        = n(1)+1       ! k is first atom in system 2
          r(:,ii)  = r(:,k)       ! Replace atom i (position ii) by atom k
          r(:,k)   = r_new        ! Put new particle in position k
          n        = n + [1,-1]   ! Move boundary up so k is now in system 1
          u(2)     = u(2) - u_old ! Update potential in system 2
          w(2)     = w(2) - w_old ! Update virial in system 2
          u(1)     = u(1) + u_new ! Update potential in system 1
          w(1)     = w(1) + w_new ! Update virial in system 1
          accepted = accepted + 1 ! Increment move counter
       end if ! End accept Metropolis test

    end if ! End choice between trying to swap 1->2 and 2->1

  end subroutine n_swap

  subroutine v_swap ( beta, dv_max, n, box, r, u, w, accepted ) ! Volume exchange move
    use maths_module,     only : metropolis, cbrt
    use potential_module, only : pot_all
    implicit none

    real(dp), intent(in)    :: beta     ! Inverse temperature
    real(dp), intent(in)    :: dv_max   ! Maximum volume change
    integer,  intent(in)    :: n(2)     ! Number of atoms in both systems
    real(dp), intent(inout) :: box(3,2) ! Box lengths in both systems
    real(dp), intent(inout) :: r(:,:)   ! Combined array of coordinates (3,n(1)+n(2))
    real(dp), intent(inout) :: u(2)     ! Total potential energy in both systems
    real(dp), intent(inout) :: w(2)     ! Total virial in both systems
    integer,  intent(inout) :: accepted ! Accepted move counter

    ! Volume exchange move in Gibbs ensemble Monte Carlo

    ! Attempts a volume exchange move, conserving the total volume
    ! Both systems are contained within the r array: r(1:n(1)) and r(n(1)+1,n(1)+n(2)).

    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values.
    ! If the move is rejected, they will be returned unchanged.

    ! NB variables such as box, n, u, w have an extra dimension, corresponding to the two simulated systems.

    real(dp) :: r_new(size(r,dim=1),size(r,dim=2)) ! New (scaled) combined coordinates

    real(dp) :: dv, delta
    real(dp) :: box_new(3,2)
    real(dp) :: v_new(2), v(2), r_scale(2), v_scale(2)
    real(dp) :: u_new(2), w_new(2)
    logical  :: overlap(2)

    if ( any ( shape(r) /= [3,sum(n)] ) ) stop 'r dimension error in v_swap'

    call random_number ( dv )                    ! Random number uniform in range (0,1)
    dv               = (2*dv-1) * dv_max         ! Now uniform in range (-dv_max,+dv_max)
    v                = product ( box, dim=1 )    ! Current volumes of both systems
    v_new            = v + [-dv,dv]              ! New volumes of both systems, total V conserved
    v_scale          = v_new / v                 ! Scaling factors for both volumes
    r_scale          = cbrt(v_scale)             ! Scaling factors for both sets of coordinates and box lengths
    box_new(:,1)     = box(:,1) * r_scale(1)     ! New box lengths, system 1
    box_new(:,2)     = box(:,2) * r_scale(2)     ! New box lengths, system 2
    r_new(:,:n(1))   = r(:,:n(1))   * r_scale(1) ! New coordinates, system 1
    r_new(:,n(1)+1:) = r(:,n(1)+1:) * r_scale(2) ! New coordinates, system 2
    call pot_all ( box_new(:,1), r_new(:,:n(1)  ), u_new(1), w_new(1), overlap(1) ) 
    call pot_all ( box_new(:,2), r_new(:,n(1)+1:), u_new(2), w_new(2), overlap(2) ) 
    if ( any ( overlap ) ) return ! Reject move on overlap

    delta = beta * sum ( u_new - u )     ! Delta U contributions for both systems
    delta = delta - n(1)*log(v_scale(1)) ! Contribution from volume scaling in system 1
    delta = delta - n(2)*log(v_scale(2)) ! Contribution from volume scaling in system 2

    if ( metropolis ( delta ) ) then ! Accept Metropolis test
       u        = u_new        ! Update potentials
       w        = w_new        ! Update virials
       box      = box_new      ! Update box lengths
       r        = r_new        ! Update coordinates
       accepted = accepted + 1 ! Increment move counter
    end if

  end subroutine v_swap

  subroutine n_move ( beta, phi, box, n, r, u, w, accepted ) ! Particle insertion or deletion move
    use maths_module,     only : random_integer, metropolis
    use potential_module, only : pot_one
    implicit none

    real(dp), intent(in)    :: beta     ! Inverse temperature
    real(dp), intent(in)    :: phi(0:)  ! Array of weights (0:nmax)
    real(dp), intent(in)    :: box(3)   ! Box lengths
    integer,  intent(inout) :: n        ! Current number of atoms
    real(dp), intent(inout) :: r(:,:)   ! Array of coordinates (3,nmax)
    real(dp), intent(inout) :: u        ! Total potential energy
    real(dp), intent(inout) :: w        ! Total virial
    integer,  intent(inout) :: accepted ! Accepted move counter

    ! Particle insertion and deletion in MultiCanonical Monte Carlo or Grand Canonical Monte Carlo
    
    ! ACKNOWLEDGEMENT
    ! An earlier code written by Leo Lue, on which this was based, is gratefully acknowledged.
    ! Needless to say, he is in no way responsible for any mistakes in this code.

    ! Carries out single-atom insertion or deletion
    ! The ensemble is defined by beta and weights phi(n)
    ! The maximum number of atoms, nmax, is defined by the dimensions of phi and r
    ! The actual number of atoms, n, must not exceed nmax

    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values.
    ! If the move is rejected, they will be returned unchanged.

    real(dp) :: r_new(3) ! Temporary position array

    integer  :: nmax, i, choice
    real(dp) :: v
    real(dp) :: u_new, w_new, u_old, w_old, delta
    logical  :: overlap

    if (   size ( r, dim=1 ) /= 3 ) stop 'r dimension error in n_move'
    nmax = size ( r, dim=2 )
    if ( ubound(phi,dim=1) /= nmax ) stop 'phi dimension error in n_move'
    if ( n > nmax  ) stop 'n error in n_move'
    v  = product ( box ) ! Volume

    choice = random_integer ( 2 ) ! Choose between insertion and deletion

    if ( choice == 1 ) then ! Attempt a particle insertion move

       if ( n == nmax ) return ! Reject: can't exceed nmax

       ! Choose random position
       call random_number ( r_new ) ! Vector with components in range (0,1)
       r_new = (2*r_new-1) * box/2  ! Components now in range (-box(:)/2,box(:)/2)

       ! Check for overlap, and calculate energy and virial, for interaction with other atoms
       call pot_one ( box, r(:,:n), r_new, 0, u_new, w_new, overlap )
       if ( overlap ) return ! Reject move on overlap

       delta = beta * u_new                         ! Delta U contribution
       delta = delta - log (v/(n+1))                ! Contribution related to (N+1)/V
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
       delta = delta + log (v/n)                    ! Contribution related to N/V
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

end module mc_module
