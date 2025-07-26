! mc_npt_module.f90
! Contains routine to carry out volume move for mc_npt program
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module mc_npt_module

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! In the "Monte Carlo at Constant Pressure" workshop, you were asked to fix
  ! the v_move routine which carries out volume moves.
  ! The calculation of delta has been FIXED in this file.

  implicit none
  private

  public :: v_move

contains

  subroutine v_move ( beta, pressure, dv_max, box, r, u, w, accepted ) ! Volume move
    use maths_module,     only : metropolis
    use potential_module, only : pot_all
    implicit none

    real,                    intent(in)    :: beta     ! Inverse temperature
    real,                    intent(in)    :: pressure ! Pressure
    real,                    intent(in)    :: dv_max   ! Maximum volume change
    real,    dimension(3),   intent(inout) :: box      ! Box lengths
    real,    dimension(:,:), intent(inout) :: r        ! Array of coordinates (3,n)
    real,                    intent(inout) :: u        ! Total potential energy
    real,                    intent(inout) :: w        ! Total virial
    integer,                 intent(inout) :: accepted ! Accepted move counter

    ! The variables declared intent(inout) contain the current (or old) values.
    ! If the move is accepted, they will be updated with the new values for the scaled system.
    ! If the move is rejected, they will be returned unchanged.

    real, dimension(size(r,dim=1),size(r,dim=2)) :: r_new ! New (scaled) coordinates

    integer            :: n
    real               :: dv, zeta, delta
    real, dimension(3) :: box_new
    real               :: v_new, v, r_scale, v_scale
    real               :: u_new, w_new
    logical            :: overlap

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in v_move'
    n  = size(r,dim=2) ! Number of atoms

    call random_number ( zeta )       ! Uniform on (0,1)
    dv      = dv_max * ( 2*zeta - 1 ) ! Delta V is uniform on (-dv_max,+dv_max)
    v       = product ( box )         ! Current (or old) volume
    v_new   = v + dv                  ! New volume
    v_scale = v_new / v               ! Scaling factor for volume
    r_scale = v_scale**(1./3.)        ! Scaling factor for coordinates and box lengths
    box_new = box * r_scale           ! New box lengths
    r_new   = r   * r_scale           ! New coordinates

    call pot_all ( box_new, r_new, u_new, w_new, overlap ) 
    if ( overlap ) return ! Reject move on overlap

    delta = beta * ( u_new - u )         ! This is the exp(-beta*Delta U) part
    delta = delta + beta * pressure * dv ! This is the exp(-beta*P*Delta V) part (corrected version)
    delta = delta - n * log ( v_scale )  ! This is the (Vnew/Vold)**N part (no correction needed)

    if ( metropolis ( delta ) ) then ! Accept Metropolis test
       u        = u_new        ! Update potential
       w        = w_new        ! Update virial
       box      = box_new      ! Update box lengths
       r        = r_new        ! Update coordinates
       accepted = accepted + 1 ! Increment move counter
    end if

  end subroutine v_move

end module mc_npt_module
