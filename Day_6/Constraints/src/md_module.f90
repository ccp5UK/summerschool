! md_module.f90
! MD algorithm & constraint routines for MD, LJ chain
module md_module

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, dp => real64

  implicit none

  ! Private variables
  private
  real(dp), parameter :: d = 0.85_dp ! Bond length for constraint routines
  real(dp), parameter :: m = 1.0_dp  ! Atomic masses (all the same)

  ! Public routines
  public :: md_intro, constraint_intro, kick, drift
  public :: rattle_a, rattle_b, milcshake_a, milcshake_b
  public :: kineng, worst_bond, thermostat

contains

  subroutine md_intro
    implicit none

    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Atomic masses', m
  
  end subroutine md_intro

  subroutine constraint_intro
    implicit none

    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Constraint bond distance', d
   
  end subroutine constraint_intro

  subroutine kick ( dt, f, p )
    implicit none
    real(dp), intent(in)    :: dt     ! Time step
    real(dp), intent(in)    :: f(:,:) ! Forces (3,n)
    real(dp), intent(inout) :: p(:,:) ! Momenta (3,n)

    ! Kick stage of velocity Verlet

    if ( any ( shape(p) /= shape(f) ) ) stop 'p dimension error in kick'

    p = p + dt * f

  end subroutine kick

  subroutine drift ( dt, p, r )
    implicit none
    real(dp), intent(in)    :: dt     ! Time step
    real(dp), intent(in)    :: p(:,:) ! Momenta (3,n)
    real(dp), intent(inout) :: r(:,:) ! Positions (3,n)

    ! Drift stage of velocity Verlet
    ! m is a module parameter declared at the top of this file

    if ( any ( shape(r) /= shape(p) ) ) stop 'r dimension error in kick'

    r = r + (dt/m) * p
    
  end subroutine drift
  
  subroutine rattle_a ( dt, q, r, p, iter )
    implicit none
    real(dp), intent(in)    :: dt     ! Time step
    real(dp), intent(in)    :: q(:,:) ! Old positions (3,n)
    real(dp), intent(inout) :: r(:,:) ! Positions (3,n)
    real(dp), intent(inout) :: p(:,:) ! Momenta (3,n)
    integer,  intent(out)   :: iter   ! Number of iterations

    ! Stage A of RATTLE, very similar to SHAKE
    ! Iteratively adjusts positions and momenta to satisfy bond constraints.

    ! On entry to this routine we assume:
    ! q stores the positions at the start of the step ("old positions")
    ! r stores the end-of-step positions following the unconstrained drift and
    ! p stores the half-step momenta following the first unconstrained half-kick
    ! All masses are the same, which slightly simplifies some of the expressions
    ! On return from this routine, r and p will hold the constrained values
    ! m and d are module parameters declared at the top of this file

    logical :: move(size(r,dim=2)), moved(size(r,dim=2)) ! Iterative flags (n)

    real(dp) :: rij(3), qij(3), dr(3)
    logical  :: done
    real(dp) :: sigma, a, lambda
    integer  :: i, j, n

    real(dp), parameter :: tol = 1.0e-9_dp,  tol2 = 2 * tol
    integer,  parameter :: iter_max = 500

    n = size(r,dim=2) ! Number of atoms
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in rattle_a'
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in rattle_a'
    if ( any ( shape(q) /= [3,n] ) ) stop 'q dimension error in rattle_a'

    iter  = 0
    done  = .false.
    moved = .true. ! Ensures that we look at each bond at least once

    do ! Iterative loop until done

       if ( done ) exit ! done is equivalent to .not.any ( moved )

       done = .true.
       move = .false.

       do i = 1, n-1 ! Loop over each constraint in turn
          j = i + 1 ! Partner atom in this constraint

          if ( moved(i) .or. moved(j) ) then ! Test whether need to re-examine ij

             rij   = r(:,i) - r(:,j)    ! Current bond vector
             sigma = d**2 - sum(rij**2) ! Amount by which constraint is violated

             if ( abs(sigma) > tol2*d**2 ) then ! Test whether constraint not already satisfied

                qij = q(:,i) - q(:,j)  ! Old bond vector

                a = 4 * dot_product ( qij, rij ) ! A

                lambda  = sigma / a        ! Linearized equation for multiplier
                dr      = qij * lambda     ! Resultant displacement vector
                p(:,i)  = p(:,i) + dr*m/dt ! Adjust i momentum
                p(:,j)  = p(:,j) - dr*m/dt ! Adjust j momentum
                r(:,i)  = r(:,i) + dr      ! Adjust i position
                r(:,j)  = r(:,j) - dr      ! Adjust j position
                move(i) = .true.           ! Flag that we moved i
                move(j) = .true.           ! Flag that we moved j
                done    = .false.          ! Flag that we moved something

             end if ! End test whether constraint not already satisfied

          end if ! End test whether need to re-examine ij

       end do ! End loop over each constraint in turn

       ! Prepare for next iteration
       moved = move
       iter  = iter + 1
       if ( iter > iter_max ) then
          write ( unit=error_unit, fmt='(a,2i10)' ) 'Too many iterations', iter, iter_max
          stop 'Error in rattle_a'
       end if

    end do ! End iterative loop until done

  end subroutine rattle_a

  subroutine rattle_b ( r, p, iter )
    implicit none
    real(dp), intent(in)    :: r(:,:) ! Positions (3,n)
    real(dp), intent(inout) :: p(:,:) ! Momenta (3,n)
    integer,  intent(out)   :: iter   ! Number of iterations

    ! Stage B of RATTLE
    ! This subroutine iteratively adjusts the momenta
    ! to satisfy the time derivatives of the bond constraints

    ! On entry to this routine we assume:
    ! r stores the positions at the end of the step with constraints applied
    ! p stores the momenta following the second unconstrained half-kick
    ! All masses are the same, which simplifies some of the expressions
    ! On return from this routine, p will hold the constrained values
    ! m and d are module parameters declared at the top of this file

    logical :: move(size(r,dim=2)), moved(size(r,dim=2)) ! Iterative flags (n)

    real(dp) :: rij(3), vij(3), delp(3)
    logical  :: done
    real(dp) :: tau, mu, b
    integer  :: i, j, n

    real(dp), parameter :: tol = 1.0e-9_dp
    integer,  parameter :: iter_max = 500

    n = size(r,dim=2) ! Number of atoms
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in rattle_b'
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in rattle_b'

    iter  = 0
    done  = .false.
    moved = .true. ! Ensures that we look at each bond at least once

    do ! Iterative loop until done

       if ( done ) exit ! done is equivalent to .not.any ( moved )

       done = .true.
       move = .false.

       do i = 1, n-1 ! Loop over each constraint in turn
          j = i + 1 ! Partner atom for this constraint

          if ( moved(i) .or. moved(j) ) then   ! Test whether need to re-examine ij
             rij = r(:,i) - r(:,j)             ! Current bond vector
             vij = p(:,i)/m - p(:,j)/m         ! Relative velocity
             tau = -m*dot_product ( rij, vij ) ! Bond derivative violation
             b   = 2*d**2                      ! B

             mu = tau / b

             if ( abs ( mu ) > tol ) then ! Test whether constraint already satisfied

                delp    = rij * mu      ! Momentum adjustment
                p(:,i)  = p(:,i) + delp ! Adjust i momentum
                p(:,j)  = p(:,j) - delp ! Adjust j momentum
                move(i) = .true.        ! Flag that we moved i
                move(j) = .true.        ! Flag that we moved j
                done    = .false.       ! Flag that we moved something

             end if ! End test whether constraint already satisfied

          end if ! End test whether need to re-examine ij

       end do ! End loop over each constraint in turn

       ! Prepare for next iteration
       moved = move
       iter  = iter + 1
       if ( iter > iter_max ) then
          write ( unit=error_unit, fmt='(a,2i10)' ) 'Too many iterations', iter, iter_max
          stop 'Error in rattle_b'
       end if

    end do ! End iterative loop until done

  end subroutine rattle_b

  subroutine milcshake_a ( dt, q, r, p, iter )
    use maths_module, only : tridiagonal
    implicit none
    real(dp), intent(in)    :: dt     ! Time step
    real(dp), intent(in)    :: q(:,:) ! Old positions (3,n)
    real(dp), intent(inout) :: r(:,:) ! Positions (3,n)
    real(dp), intent(inout) :: p(:,:) ! Momenta (3,n)
    integer,  intent(out)   :: iter   ! Number of iterations

    ! This subroutine iteratively adjusts the positions stored in the array r
    ! and the momenta stored in the array p, to satisfy the bond constraints
    ! using a tri-diagonal solver.
    ! See AG Bailey, CP Lowe, and AP Sutton, J Comput Phys, 227, 8949 (2008)
    ! and AG Bailey, CP Lowe, and AP Sutton, Comput Phys Commun, 180, 594 (2009)

    ! On entry to this routine we assume:
    ! q stores the positions at the start of the step
    ! r stores the positions following the unconstrained drift and
    ! p stores the momenta following the first unconstrained half-kick
    ! All masses are the same, which simplifies some of the expressions
    ! On return from this routine, r and p will hold the constrained values
    ! and iter will store the number of iterations carried out to converge
    ! m and d are module parameters declared at the top of this file

    integer :: n ! number of atoms
    integer :: c ! number of constraints, n-1

    real(dp) :: qij(3,size(r,dim=2)-1) ! old bond vectors (3,c)
    real(dp) :: rij(3,size(r,dim=2)-1) ! new bond vectors (3,c)
    real(dp) :: dr(3,size(r,dim=2)-1)  ! position update vectors (3,c)

    real(dp) :: lambda(size(r,dim=2)-1) ! multipliers (c)
    real(dp) :: sigma(size(r,dim=2)-1)  ! amounts by which constraints are violated (c)
    real(dp) :: ad(size(r,dim=2)-1)     ! diagonal elements of A matrix (c)
    real(dp) :: au(size(r,dim=2)-2)     ! upper-diagonal elements of A matrix (c-1)
    real(dp) :: al(size(r,dim=2)-2)     ! lower-diagonal elements of A matrix (c-1)

    real(dp), parameter :: tol = 1.0e-9_dp, tol2 = 2.0_dp * tol
    integer,  parameter :: iter_max = 500

    n = size(r,dim=2) ! Number of atoms
    c = n - 1         ! Number of constraints
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in milcshake_a'
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in milcshake_a'
    if ( any ( shape(q) /= [3,n] ) ) stop 'q dimension error in milcshake_a'

    qij = q(:,1:c) - q(:,2:n) ! Old bond vectors
    rij = r(:,1:c) - r(:,2:n) ! New (unconstrained) bond vectors

    ! Elements of tridiagonal matrix A (dot products of old and new bond vectors)
    al(1:c-1) = -2*sum( qij(:,1:c-1)*rij(:,2:c),   dim=1 ) ! c-1 elements of lower-diagonal
    ad(1:c)   =  4*sum( qij(:,1:c)  *rij(:,1:c),   dim=1 ) ! c elements of diagonal
    au(1:c-1) = -2*sum( qij(:,2:c)  *rij(:,1:c-1), dim=1 ) ! c-1 elements of upper-diagonal

    iter = 0

    do ! Iterative loop until done

       rij   = r(:,1:c) - r(:,2:n)      ! New bond vectors
       sigma = d**2 - sum(rij**2,dim=1) ! Amounts by which constraints are violated

       if ( maxval(abs(sigma)) <= tol2*d**2 ) exit ! Test for done

       lambda = tridiagonal ( al, ad, au, sigma )  ! Solve tridiagonal equation
       dr     = spread(lambda,dim=1,ncopies=3)*qij ! Resultant update vectors

       p(:,1:c) = p(:,1:c) + dr*m/dt ! Adjust i momenta
       p(:,2:n) = p(:,2:n) - dr*m/dt ! Adjust j momenta
       r(:,1:c) = r(:,1:c) + dr      ! Adjust i positions
       r(:,2:n) = r(:,2:n) - dr      ! Adjust j positions

       iter = iter + 1
       if ( iter > iter_max ) then
          write ( unit=error_unit, fmt='(a,2i15)' ) 'Too many iterations', iter, iter_max
          stop 'Error in milcshake_a'
       end if

    end do ! End iterative loop until done

  end subroutine milcshake_a

  subroutine milcshake_b ( r, p, iter )
    use maths_module, only : tridiagonal
    implicit none
    real(dp), intent(in)    :: r(:,:) ! Positions (3,n)
    real(dp), intent(inout) :: p(:,:) ! Momenta (3,n)
    integer,  intent(out)   :: iter   ! Number of iterations (always 1)

    ! This subroutine adjusts the momenta stored in the array p
    ! to satisfy the time derivatives of the bond constraints using a tridiagonal solver.
    ! See AG Bailey, CP Lowe, and AP Sutton, J Comput Phys, 227, 8949 (2008)
    ! and AG Bailey, CP Lowe, and AP Sutton, Comput Phys Commun, 180, 594 (2009)

    ! On entry to this routine we assume:
    ! r stores the positions at the end of the step with constraints applied
    ! p stores the momenta following the second unconstrained half-kick
    ! All masses are the same, which simplifies some of the expressions
    ! On return from this routine, p will hold the constrained values
    ! m is a module parameter declared at the top of this file

    integer :: n ! number of atoms
    integer :: c ! number of constraints, n-1

    real(dp) :: rij(3,size(r,dim=2)-1)  ! bond vectors (3,c)
    real(dp) :: vij(3,size(r,dim=2)-1)  ! new relative velocities (3,c)
    real(dp) :: delp(3,size(r,dim=2)-1) ! momentum change vectors (3,c)

    real(dp) :: mu(size(r,dim=2)-1)  ! multipliers (c)
    real(dp) :: tau(size(r,dim=2)-1) ! amounts by which constraints are violated (c)
    real(dp) :: bd(size(r,dim=2)-1)  ! diagonal elements of B matrix (c)
    real(dp) :: bu(size(r,dim=2)-2)  ! upper-diagonal elements of B matrix (c-1)
    real(dp) :: bl(size(r,dim=2)-2)  ! lower-diagonal elements of B matrix (c-1)

    n = size(r,dim=2) ! Number of atoms
    c = n - 1         ! Number of constraints
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in milcshake_b'
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in milcshake_b'

    rij = r(:,1:c) - r(:,2:n) ! Bond vectors

    ! Elements of tridiagonal matrix B (dot products of bond vectors)
    bl =  -sum ( rij(:,1:c-1)*rij(:,2:c),   dim=1 ) ! c-1 elements of lower-diagonal
    bd = 2*sum ( rij(:,1:c)  *rij(:,1:c),   dim=1 ) ! c elements of diagonal
    bu =  -sum ( rij(:,2:c)  *rij(:,1:c-1), dim=1 ) ! c-1 elements of upper-diagonal

    vij  = p(:,1:c)/m - p(:,2:n)/m         ! Relative velocities
    tau  = -m*sum( vij*rij, dim=1 )        ! Amounts by which constraints are violated
    mu   = tridiagonal ( bl, bd, bu, tau ) ! Solve tridiagonal system
    delp = spread(mu,dim=1,ncopies=3)*rij  ! Resultant update vectors

    p(:,1:c) = p(:,1:c) + delp ! Adjust i momenta
    p(:,2:n) = p(:,2:n) - delp ! Adjust j momenta

    iter = 1

  end subroutine milcshake_b

  subroutine thermostat ( temperature, nfree, p )
    use maths_module, only : metropolis
    implicit none
    real(dp), intent(in)    :: temperature ! Temperature
    integer,  intent(in)    :: nfree       ! Number of degrees of freedom
    real(dp), intent(inout) :: p(:,:)      ! Momenta (3,n)

    ! Implements the thermostat of DM Heyes, Chem. Phys., 82, 285 (1983)
    ! corrected by M Hecht et al, Phys. Rev. E, 72, 011408 (2005),
    ! using the method suggested by D Frenkel, B Smit, Understanding Molecular Simulation, 3rd ed (2023),
    ! with an erratum at https://github.com/UnderstandingMolecularSimulation/List-of-Errata-and-Additions

    ! All masses are the same
    ! m is a module parameter declared at the top of this file

    ! Notes:
    ! (1) Rescaling momenta this way only conserves variables such as total momentum and relative velocity 
    !     components along bond constraints, in the sense that they remain =0 if already =0.
    ! (2) nfree must contain the number of degrees of freedom equal to the
    !     number of momentum components 3N, minus the number of conserved variables and constraints.
    ! (3) We have set, arbitrarily, the maximum change in lnz, del_lnz, to a fairly small value,
    !     but more generally it would be a simulation parameter, adjustable by the user.

    real(dp) :: lnz, z, delta
    integer  :: n

    real(dp), parameter :: del_lnz = 0.01_dp 

    n =  size(p,dim=2)
    if ( size(p,dim=1) /= 3 ) stop 'p dimension error in andersen'

    call random_number(lnz)
    lnz = (2.0_dp*lnz - 1.0_dp) * del_lnz ! Choose ln(z) in range +/- del_lnz
    z   = exp(lnz)                        ! Momentum scaling factor

    delta = -nfree*lnz + kineng(p)*(z**2-1.0_dp)/temperature
    if ( metropolis(delta) ) then
       p = p * z
    end if

  end subroutine thermostat

  function kineng ( p ) result (k)
    implicit none
    real(dp), intent(in) :: p(:,:) ! Momenta (3,n)
    real(dp)             :: k      ! Total kinetic energy

    ! Returns total kinetic energy
    ! All masses are the same
    ! m is a module parameter declared at the top of this file

    k = sum(p**2)/(2*m)

  end function kineng

  subroutine worst_bond ( r, p )
    use maths_module, only : pbc
    implicit none
    real(dp), intent(in) :: r(:,:) ! Positions (3,n)
    real(dp), intent(in) :: p(:,:) ! Momenta (3,n)

    ! Writes out 2 measures relating to bond constraints between neighbouring atoms
    ! The "deviation" is |dij-d| where dij is the magnitude of the vector rij = ri-rj
    ! The "derivative" is the time derivative of this, which turns out to be |vij.rij|/dij
    ! where vij=vi-vj is the relative velocity vector.
    ! In each case, we print the worst values for the whole molecule.
    ! Both these quantities would be zero for perfectly constrained bonds,
    ! but neither of them is dimensionless, and they are just roughly indicative of what is (nearly) conserved.
    ! All masses are the same
    ! m and d are module parameters declared at the top of this file

    integer :: n ! Number of atoms
    integer :: c ! Number of constraints, n-1

    real(dp) :: rij(3,size(r,dim=2)-1) ! bond vectors (3,c)
    real(dp) :: vij(3,size(r,dim=2)-1) ! relative velocities (3,c)

    real(dp) :: dij(size(r,dim=2)-1)   ! bond lengths (c)
    real(dp) :: sigma(size(r,dim=2)-1) ! bond length deviations (c)
    real(dp) :: tau(size(r,dim=2)-1)   ! bond length derivatives (c)

    n = size(r,dim=2) ! Number of atoms
    c = n - 1         ! Number of constraints
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in worst_bond'
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in worst_bond'

    rij   = r(:,1:c) - r(:,2:n)           ! Relative position vectors
    dij   = norm2(rij,dim=1)              ! Distances between atoms
    sigma = abs(dij-d)                    ! Deviations from specified bond lengths
    vij   = p(:,1:c)/m - p(:,2:n)/m       ! Relative velocity vectors
    tau   = abs(sum(rij*vij,dim=1)) / dij ! Time derivatives

    write ( unit=output_unit, fmt='(a,t31,es14.4)') 'Worst bond length deviation ', maxval(sigma)
    write ( unit=output_unit, fmt='(a,t31,es14.4)') 'Worst bond length derivative', maxval(tau)

  end subroutine worst_bond

end module md_module
