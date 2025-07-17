! mc_gibbs_module.f90
! Routines to carry out particle exchange and volume exchange moves for mc_gibbs program
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module mc_gibbs_module

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  implicit none
  private

  public :: n_move, v_move

contains

  subroutine n_move ( beta, box, n, r, u, w, accepted ) ! Particle exchange move
    use maths_module,     only : random_integer, metropolis
    use potential_module, only : pot_one
    implicit none

    real,                    intent(in)    :: beta     ! Inverse temperature
    real,    dimension(3,2), intent(in)    :: box      ! Box lengths in both systems
    integer, dimension(2),   intent(inout) :: n        ! Number of atoms in both systems
    real,    dimension(:,:), intent(inout) :: r        ! Array of coordinates (3,n(1)+n(2))
    real,    dimension(2),   intent(inout) :: u        ! Total potential energy in both systems
    real,    dimension(2),   intent(inout) :: w        ! Total virial in both systems
    integer,                 intent(inout) :: accepted ! Accepted move counter

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

    integer              :: i, k, ii, direction
    real,   dimension(3) :: r_new
    real,   dimension(2) :: v
    real                 :: u_new, w_new, u_old, w_old, delta
    logical              :: overlap

    if ( any ( shape(r) /= [3,sum(n)] ) ) stop 'r dimension error in n_move'

    v = product ( box, dim=1 ) ! Both box volumes

    direction = random_integer(2) ! Choose direction of swap

    if ( direction == 2 ) then ! Try swapping 1->2

       if ( n(1) <= 1 ) return ! Disallow n(1)->0

       i = random_integer ( n(1) ) ! Choose atom in system 1 at random
       call pot_one ( box(:,1), r(:,:n(1)), r(:,i), i, u_old, w_old, overlap )
       if ( overlap ) stop 'Overlap in current configuration' ! Should never happen
       call random_number ( r_new )   ! 3 random numbers in range (0,1)
       r_new = (r_new - 0.5)*box(:,2) ! Random position in system 2
       call pot_one ( box(:,2), r(:,n(1)+1:), r_new, 0, u_new, w_new, overlap )
       if ( overlap ) return ! Reject move on overlap

       delta = beta * ( u_new - u_old )               ! Delta U contribution
       delta = delta - log ( v(2) / real ( n(2)+1 ) ) ! Contribution from creation in 2
       delta = delta + log ( v(1) / real ( n(1)   ) ) ! Contribution from destruction in 1

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
       call random_number ( r_new )   ! 3 random numbers in range (0,1)
       r_new = (r_new - 0.5)*box(:,1) ! Random position in system 1
       call pot_one ( box(:,1), r(:,:n(1)), r_new, 0, u_new, w_new, overlap )
       if ( overlap ) return ! Reject move on overlap

       delta = beta * ( u_new - u_old )               ! Delta U contribution
       delta = delta - log ( v(1) / real ( n(1)+1 ) ) ! Contribution from creation in 1
       delta = delta + log ( v(2) / real ( n(2)   ) ) ! Contribution from destruction in 2

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

  end subroutine n_move

  subroutine v_move ( beta, dv_max, n, box, r, u, w, accepted ) ! Volume exchange move
    use maths_module,     only : metropolis
    use potential_module, only : pot_all
    implicit none

    real,                    intent(in)    :: beta     ! Inverse temperature
    real,                    intent(in)    :: dv_max   ! Maximum volume change
    integer, dimension(2),   intent(in)    :: n        ! Number of atoms in both systems
    real,    dimension(3,2), intent(inout) :: box      ! Box lengths in both systems
    real,    dimension(:,:), intent(inout) :: r        ! Array of coordinates (3,n(1)+n(2))
    real,    dimension(2),   intent(inout) :: u        ! Total potential energy in both systems
    real,    dimension(2),   intent(inout) :: w        ! Total virial in both systems
    integer,                 intent(inout) :: accepted ! Accepted move counter

    ! Attempts a volume exchange move, conserving the total volume
    ! Both systems are contained within the r array: r(1:n(1)) and r(n(1)+1,n(1)+n(2)).
    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values.
    ! If the move is rejected, they will be returned unchanged.

    real, dimension(size(r,dim=1),size(r,dim=2)) :: r_new ! New (scaled) coordinates in both systems

    real                    :: dv, zeta, delta
    real,    dimension(3,2) :: box_new
    real,    dimension(2)   :: v_new, v, r_scale, v_scale
    real,    dimension(2)   :: u_new, w_new
    logical, dimension(2)   :: overlap

    if ( any ( shape(r) /= [3,sum(n)] ) ) stop 'r dimension error in v_move'

    call random_number ( zeta )                  ! Uniform on (0,1)
    dv               = dv_max * ( 2*zeta - 1 )   ! Delta V uniform on (-dv_max,+dv_max)
    v                = product ( box, dim=1 )    ! Current volumes of both systems
    v_new            = v + [-dv,dv]              ! New volumes of both systems, total V conserved
    v_scale          = v_new / v                 ! Scaling factors for both volumes
    r_scale          = v_scale**(1./3.)          ! Scaling factors for coordinates and box lengths
    box_new(:,1)     = box(:,1) * r_scale(1)     ! New box lengths, system 1
    box_new(:,2)     = box(:,2) * r_scale(2)     ! New box lengths, system 2
    r_new(:,:n(1))   = r(:,:n(1))   * r_scale(1) ! New coordinates, system 1
    r_new(:,n(1)+1:) = r(:,n(1)+1:) * r_scale(2) ! New coordinates, system 2
    call pot_all ( box_new(:,1), r_new(:,:n(1)  ), u_new(1), w_new(1), overlap(1) ) 
    call pot_all ( box_new(:,2), r_new(:,n(1)+1:), u_new(2), w_new(2), overlap(2) ) 
    if ( any ( overlap ) ) return ! Reject move on overlap

    delta = beta * sum ( u_new - u )           ! Delta U contributions for both systems
    delta = delta - real(n(1))*log(v_scale(1)) ! Contribution from volume scaling in system 1
    delta = delta - real(n(2))*log(v_scale(2)) ! Contribution from volume scaling in system 2

    if ( metropolis ( delta ) ) then ! Accept Metropolis test
       u        = u_new        ! Update potentials
       w        = w_new        ! Update virials
       box      = box_new      ! Update box lengths
       r        = r_new        ! Update coordinates
       accepted = accepted + 1 ! Increment move counter
    end if

  end subroutine v_move

end module mc_gibbs_module
