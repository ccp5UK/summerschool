! force_module.f90
! Force routines for MD, WCA-LJ chain
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module force_module

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit

  implicit none
  private

  ! Private variables
  real, parameter :: r_cut    = 2.0**(1.0/6.0) ! WCA LJ potential cutoff
  real, parameter :: r_cut_sq = r_cut**2       ! and its square
  real, parameter :: d        = 2.0**(1.0/6.0) ! Bond length for spring potential
  real, parameter :: kappa    = 10000.0        ! Force constant for spring potential

  ! Public routines
  public :: force_intro, force, spring

contains

  subroutine force_intro ( springs )
    implicit none
    logical, intent(in) :: springs ! .true. for model with springs, .false. for constraints

    write ( unit=output_unit, fmt='(a)'           ) 'WCA Lennard-Jones chain'   
    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Diameter sigma',     1.0
    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Well depth epsilon', 1.0
    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Potential cutoff',   r_cut
    if ( springs ) then
       write ( unit=output_unit, fmt='(a,t31,f10.4)'  ) 'Spring bond distance', d
       write ( unit=output_unit, fmt='(a,t31,es14.4)' ) 'Force constant kappa', kappa
    end if
   
  end subroutine force_intro

  subroutine force ( box, r, u, f )
    use maths_module, only : pbc
    implicit none
    real, dimension(3),   intent(in)  :: box ! Box lengths
    real, dimension(:,:), intent(in)  :: r   ! Positions (3,n)
    real,                 intent(out) :: u   ! Potential energy
    real, dimension(:,:), intent(out) :: f   ! Forces (3,n)

    ! u is the WCA LJ potential energy for the atomic chain (omitting bonded neighbours)
    ! The Lennard-Jones energy and sigma parameters are taken to be epsilon = 1, sigma = 1
    ! Positions are assumed to be in these units
    ! Forces are calculated in the same units and stored in the array f
    ! r_cut and r_cut_sq are module parameters declared at the top of this file.

    integer              :: i, j, n
    real                 :: rij_sq, sr2, sr6, sr12, uij
    real, dimension(3)   :: rij, fij
    real, parameter      :: sr2_ovr = 1.77 ! Overlap threshold (uij > 100)

    n = size(r,dim=2) ! Number of atoms
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in force'
    if ( any ( shape(f) /= [3,n] ) ) stop 'f dimension error in force'
    if ( minval(box) < 2.0*r_cut   ) stop 'Box too small'

    ! Initialize
    f = 0.0
    u = 0.0

    do i = 1, n - 2 ! Begin outer loop over atoms, stopping 2 short of the end

       do j = i + 2, n ! Begin inner loop over atoms omitting nearest neighbour

          rij    = r(:,i) - r(:,j)       ! Separation vector
          rij    = pbc ( rij, box )      ! Periodic boundary correction
          rij_sq = sum ( rij**2 )        ! Squared separation
          if ( rij_sq > r_cut_sq ) cycle ! Skip if outside cutoff

          sr2 = 1.0 / rij_sq  ! (sigma/rij)**2
          if ( sr2 > sr2_ovr ) stop 'overlap detected'

          sr6  = sr2 ** 3
          sr12 = sr6 ** 2
          uij  =  4 * (   sr12 - sr6 ) + 1         ! WCA LJ pair potential
          fij  = 24 * ( 2*sr12 - sr6 ) * sr2 * rij ! LJ Pair forces

          u      = u      + uij
          f(:,i) = f(:,i) + fij
          f(:,j) = f(:,j) - fij

       end do ! End inner loop over atoms

    end do ! End outer loop over atoms

  end subroutine force

  subroutine spring ( box, r, v, g )
    use maths_module, only : pbc
    implicit none
    real, dimension(3),   intent(in)  :: box ! Box lengths
    real, dimension(:,:), intent(in)  :: r   ! Positions (3,n)
    real,                 intent(out) :: v   ! Total harmonic spring potential energy
    real, dimension(:,:), intent(out) :: g   ! Forces (3,n)

    ! Calculates bond spring potential for atomic chain
    ! Forces are also calculated and stored in the array g
    ! d and kappa are module parameters declared at the top of this file.

    integer            :: i, j, n
    real               :: dij, vij
    real, dimension(3) :: rij, gij

    n = size(r,dim=2) ! Number of atoms
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in spring'
    if ( any ( shape(g) /= [3,n] ) ) stop 'g dimension error in spring'

    ! Initialize
    g = 0.0
    v = 0.0

    do i = 1, n - 1 ! Begin loop over bonds
       j = i + 1    ! Nearest neighbour

       rij = r(:,i) - r(:,j)  ! Separation vector
       rij = pbc ( rij, box ) ! Periodic boundary correction
       dij = norm2(rij)       ! Separation

       vij = 0.5 * kappa * (dij-d) ** 2      ! Spring pair potential
       gij = kappa * rij * ( d - dij ) / dij ! Spring pair force

       v      = v      + vij
       g(:,i) = g(:,i) + gij
       g(:,j) = g(:,j) - gij

    end do ! End loop over bonds

  end subroutine spring

end module force_module
