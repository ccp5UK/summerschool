! force_module.f90
! Force routines for MD, LJ chain
module force_module

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, dp => real64

  implicit none

  ! Private variables
  private
  real(dp), parameter :: d        = 0.85_dp    ! Bond length for spring potential
  real(dp), parameter :: kappa    = 10000.0_dp ! Force constant for spring potential

  ! Public routines
  public :: force_intro, spring_intro, force, spring

contains

  subroutine force_intro
    implicit none

    write ( unit=output_unit, fmt='(a)'           ) 'Lennard-Jones chain'   
    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Diameter sigma',     1.0_dp
    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Well depth epsilon', 1.0_dp
    write ( unit=output_unit, fmt='(a,t31)'       ) 'NO potential cutoff'

  end subroutine force_intro

  subroutine spring_intro
    implicit none

    write ( unit=output_unit, fmt='(a,t31,f10.4)'  ) 'Spring bond distance', d
    write ( unit=output_unit, fmt='(a,t31,es14.4)' ) 'Force constant kappa', kappa

  end subroutine spring_intro

  subroutine force ( r, u, f )
    implicit none
    real(dp), intent(in)  :: r(:,:) ! Positions (3,n)
    real(dp), intent(out) :: u      ! Potential energy
    real(dp), intent(out) :: f(:,:) ! Forces (3,n)

    ! u is the LJ potential energy for the atomic chain (omitting bonded neighbours)
    ! The Lennard-Jones energy and sigma parameters are taken to be epsilon = 1, sigma = 1
    ! Positions are assumed to be in these units
    ! Forces are calculated in the same units and stored in the array f
    ! There is no periodic box and no potential cutoff

    integer  :: i, j, n
    real(dp) :: rij_sq, sr2, sr6, sr12, uij
    real(dp) :: rij(3), fij(3)

    real(dp), parameter :: sr2_ovr = 1.77_dp ! Overlap threshold (uij > 100)

    n = size(r,dim=2) ! Number of atoms
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in force'
    if ( any ( shape(f) /= [3,n] ) ) stop 'f dimension error in force'

    ! Initialize
    f = 0
    u = 0

    do i = 1, n - 2 ! Begin outer loop over atoms, stopping 2 short of the end

       do j = i + 2, n ! Begin inner loop over atoms omitting nearest neighbour

          rij    = r(:,i) - r(:,j)       ! Separation vector
          rij_sq = sum ( rij**2 )        ! Squared separation

          sr2 = 1 / rij_sq  ! (sigma/rij)**2
          if ( sr2 > sr2_ovr ) stop 'overlap detected'

          sr6  = sr2 ** 3
          sr12 = sr6 ** 2
          uij  =  4 * (   sr12 - sr6 )             ! LJ pair potential
          fij  = 24 * ( 2*sr12 - sr6 ) * sr2 * rij ! LJ Pair forces

          u      = u      + uij
          f(:,i) = f(:,i) + fij
          f(:,j) = f(:,j) - fij

       end do ! End inner loop over atoms

    end do ! End outer loop over atoms

  end subroutine force

  subroutine spring ( r, v, g )
    implicit none
    real(dp), intent(in)  :: r(:,:) ! Positions (3,n)
    real(dp), intent(out) :: v      ! Total harmonic spring potential energy
    real(dp), intent(out) :: g(:,:) ! Forces (3,n)

    ! Calculates bond spring potential for atomic chain
    ! Forces are also calculated and stored in the array g
    ! d and kappa are module parameters declared at the top of this file.

    integer  :: i, j, n
    real(dp) :: dij, vij
    real(dp) :: rij(3), gij(3)

    n = size(r,dim=2) ! Number of atoms
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in spring'
    if ( any ( shape(g) /= [3,n] ) ) stop 'g dimension error in spring'

    ! Initialize
    g = 0
    v = 0

    do i = 1, n - 1 ! Begin loop over bonds
       j = i + 1    ! Nearest neighbour

       rij = r(:,i) - r(:,j)  ! Separation vector
       dij = norm2(rij)       ! Separation

       vij = (kappa/2) * (dij-d) ** 2        ! Spring pair potential
       gij = kappa * rij * ( d - dij ) / dij ! Spring pair force

       v      = v      + vij
       g(:,i) = g(:,i) + gij
       g(:,j) = g(:,j) - gij

    end do ! End loop over bonds

  end subroutine spring

end module force_module
