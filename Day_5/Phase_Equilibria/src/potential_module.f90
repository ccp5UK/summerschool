! potential_module.f90
! Potential energy routines for cut and shifted LJ potential
module potential_module

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, dp => real64

  implicit none

  ! Private variables
  private
  real(dp), parameter :: r_cut = 2.5_dp                                  ! Potential cutoff
  real(dp), parameter :: r_cut_sq = r_cut**2                             ! and its square
  real(dp), parameter :: u_cut = 4 * ( r_cut_sq**(-6) - r_cut_sq**(-3) ) ! Potential at cutoff

  ! Public routines
  public :: potential_intro, pot_one, pot_all, t_config

contains

  subroutine potential_intro
    implicit none

    write ( unit=output_unit, fmt='(a)'           ) 'Lennard-Jones cut-and-shifted potential'
    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Diameter sigma',     1.0_dp
    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Well depth epsilon', 1.0_dp
    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Potential cutoff',   r_cut

  end subroutine potential_intro
  
  subroutine pot_all ( box, r, u, w, overlap )
    use maths_module, only : pbc
    implicit none

    real(dp), intent(in)  :: box(3)  ! Simulation box lengths
    real(dp), intent(in)  :: r(:,:)  ! Atom coordinates (3,n)
    real(dp), intent(out) :: u       ! Total potential energy
    real(dp), intent(out) :: w       ! Total virial
    logical,  intent(out) :: overlap ! Flag for overlap

    ! Calculates sum of potential interactions between all pairs in the system
    ! Returns total potential u and virial w, unless an overlap is detected
    ! in which case this is flagged with overlap=.true., and u, w, should not be used
    ! r_cut and r_cut_sq are module parameters declared at the top of this file
    
    integer  :: i, j, n
    real(dp) :: u_ij, w_ij, rij_sq
    real(dp) :: rij(3)

    if ( size ( r, dim=1 ) /= 3 ) stop 'r dimension error in pot_all'
    n  = size ( r, dim=2 ) ! Number of atoms

    ! Check all box dimensions against potential cutoff
    ! This test is typically made at the start of a run and also following a volume scaling attempt
    ! In either case, it is best to stop the program, and check things
    if ( minval(box) < 2*r_cut ) then
       write ( unit=error_unit, fmt='(a,3f10.4)') 'Box too small in pot_all', box
       stop 'Error in pot_all'
    end if

    ! Initialize
    u = 0
    w = 0
    overlap = .false.

    ! Quick exit if n is too small
    if ( n < 2 ) return
    
    ! Double loop over all distinct pairs
    do i = 1, n-1
       do j = i+1, n

          rij    = r(:,i) - r(:,j)  ! Separation vector, to which we apply
          rij    = pbc ( rij, box ) ! the minimum image correction
          rij_sq = sum ( rij**2 )   ! Squared separation

          if ( rij_sq < r_cut_sq ) then ! Check within range

             call pot_pair ( rij_sq, u_ij, w_ij, overlap )
             if ( overlap ) return

             u = u + u_ij ! Accumulate potential energy
             w = w + w_ij ! and total virial

          end if ! End check within range

       end do
    end do
    ! End double loop over all distinct pairs

  end subroutine pot_all

  subroutine pot_one ( box, r, r_i, i, u, w, overlap )
    use maths_module, only : pbc
    implicit none

    real(dp), intent(in)  :: box(3)  ! Simulation box lengths
    real(dp), intent(in)  :: r(:,:)  ! Atom coordinates (3,n)
    real(dp), intent(in)  :: r_i(3)  ! Coordinates of atom of interest
    integer,  intent(in)  :: i       ! Index to be skipped
    real(dp), intent(out) :: u       ! Potential energy
    real(dp), intent(out) :: w       ! Virial
    logical,  intent(out) :: overlap ! Flag for overlap

    ! Calculates sum of potential interactions between one atom, r_i, and all others.
    ! In the case that r_i represents the coordinates of atom i, we want to skip j=i.
    ! However, r_i may not correspond to an existing atom, in which case i=0.
    ! Returns potential u and virial w for this atom, unless overlap is detected
    ! in which case this is flagged with overlap=.true., and u, w, should not be used
    ! r_cut and r_cut_sq are module parameters declared at the top of this file
 
    integer   :: j, n
    real(dp)  :: rij_sq, u_ij, w_ij
    real(dp)  :: rij(3)

    if ( size ( r, dim=1 ) /= 3 ) stop 'r dimension error in pot_one'
    n  = size ( r, dim=2 ) ! Number of atoms

    ! Check all box dimensions against potential cutoff
    if ( minval(box) < 2*r_cut ) then ! This should never happen
       write ( unit=error_unit, fmt='(a,3f10.4)') 'Box too small in pot_one', box
       stop 'Error in pot_one'
    end if

    ! Initialize
    u = 0
    w = 0
    overlap = .false.

    ! Quick exit if n is too small
    if ( n < 1 ) return

    do j = 1, n ! Loop over all atoms

       if ( i == j ) cycle ! Skip self

       rij    = r_i - r(:,j)     ! Separation vector, to which we apply
       rij    = pbc ( rij, box ) ! the minimum image correction
       rij_sq = sum ( rij**2 )   ! Squared separation

       if ( rij_sq < r_cut_sq ) then ! Check within range

          call pot_pair ( rij_sq, u_ij, w_ij, overlap )
          if ( overlap ) return

          u = u + u_ij ! Accumulate potential energy
          w = w + w_ij ! and virial

       end if ! End check within range

    end do ! End loop over all atoms

  end subroutine pot_one

  subroutine pot_pair ( rsq, u, w, overlap )
    implicit none

    real(dp), intent(in)  :: rsq     ! Squared separation
    real(dp), intent(out) :: u       ! Pair potential, shifted at cutoff
    real(dp), intent(out) :: w       ! Pair virial
    logical,  intent(out) :: overlap ! Flag for overlap

    real(dp)            :: sr2, sr6, sr12
    real(dp), parameter :: sr2_max = 1.77_dp ! overlap threshold (u > 100)

    ! Computes pair potential and virial for given squared separation
    ! If overlap is detected it is flagged with overlap=.true.
    ! It is assumed that rsq was calculated taking periodic boundaries into account
    ! and that this routine is only called if rsq < r_cut_sq
    ! This is the shifted Lennard-Jones potential, with u_cut as a module parameter
    ! r_cut_sq and u_cut are module parameters declared at the top of this file

    sr2      = 1 / rsq
    overlap  = sr2 > sr2_max

    sr6  = sr2**3
    sr12 = sr6**2
    u    = 4  * ( sr12 - sr6 ) - u_cut
    w    = 24 * ( 2*sr12 - sr6  )

  end subroutine pot_pair

  function t_config ( box, r ) result ( t )
    use maths_module, only : pbc
    implicit none

    real(dp), intent(in)  :: box(3) ! Simulation box lengths
    real(dp), intent(in)  :: r(:,:) ! Atom coordinates (3,n)
    real(dp)              :: t      ! Configurational temperature estimate

    ! Configurational temperature calculation
    ! This is the squared force on each atom
    ! divided by the Laplacian of the potential for each atom
    ! We sum the numerator and denominator over atoms to give a single result, t.
    ! r_cut and r_cut_sq are module parameters declared at the top of this file

    real(dp) :: rij(3), fij(3)
    real(dp) :: f(size(r,dim=1),size(r,dim=2))
    
    integer  :: i, j, n
    real(dp) :: rij_sq, sr2, sr6, sr12, lap

    if ( size ( r, dim=1 ) /= 3 ) stop 'r dimension error in t_config'
    n  = size ( r, dim=2 ) ! Number of atoms

    ! We require that all box dimensions are sufficiently large
    if ( minval(box) < 2*r_cut ) then ! This should never happen
       write ( unit=error_unit, fmt='(a,3f10.4)') 'Box too small in t_config', box(:)
       stop 'Error in t_config'
    end if

    ! Initialize
    lap = 0
    f   = 0

    ! Double loop over all distinct pairs
    do i = 1, n-1
       do j = i+1, n

          rij    = r(:,i) - r(:,j)   ! Separation vector
          rij    = pbc ( rij, box )  ! Periodic boundaries
          rij_sq = sum ( rij**2 )    ! Squared separation

          if ( rij_sq < r_cut_sq ) then ! Check within range

             sr2    = 1 / rij_sq
             sr6    = sr2**3
             sr12   = sr6**2
             fij    = rij * (2*sr12 - sr6) *sr2       ! LJ pair forces
             f(:,i) = f(:,i) + fij                    ! Increment force on i
             f(:,j) = f(:,j) - fij                    ! Increment force on j
             lap    = lap + ( 22*sr12 - 5*sr6 ) * sr2 ! LJ pair Laplacian

          end if ! End check within range

       end do
    end do
    ! End double loop over all distinct pairs

    lap = lap * 24 * 2       ! Numerical factor 24*epsilon and factor 2 for ij and ji
    f   = f * 24             ! Numerical factor 24*epsilon
    t   = sum ( f**2 ) / lap ! Result

  end function t_config

end module potential_module
