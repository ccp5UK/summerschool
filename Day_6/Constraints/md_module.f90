! md_module.f90
! Force & constraint routines for MD, WCA-LJ chain
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module md_module

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit

  implicit none
  private

  ! Private variables
  real, parameter :: d = 2.0**(1.0/6.0) ! Bond length for constraint routines
  real, parameter :: m = 1              ! Atomic masses (all the same)

  ! Public routines
  public :: md_intro, kick, drift
  public :: rattle_a, rattle_b, milcshake_a, milcshake_b
  public :: kineng, worst_bond, andersen

contains

  subroutine md_intro ( springs )
    implicit none
    logical, intent(in) :: springs ! .true. for model with springs, .false. for constraints

    write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Atomic masses', m
    if ( .not. springs ) then
       write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Constraint bond distance', d
    end if
   
  end subroutine md_intro

  subroutine kick ( dt, f, p )
    implicit none
    real,                 intent(in)    :: dt ! Time step
    real, dimension(:,:), intent(in)    :: f  ! Forces (3,n)
    real, dimension(:,:), intent(inout) :: p  ! Momenta (3,n)

    ! Kick stage of velocity Verlet

    if ( any ( shape(p) /= shape(f) ) ) stop 'p dimension error in kick'

    p = p + dt * f

  end subroutine kick

  subroutine drift ( dt, p, r )
    implicit none
    real,                 intent(in)    :: dt ! Time step
    real, dimension(:,:), intent(in)    :: p  ! Momenta (3,n)
    real, dimension(:,:), intent(inout) :: r  ! Positions (3,n)

    ! Drift stage of velocity Verlet
    ! m is a module parameter declared at the top of this file

    if ( any ( shape(r) /= shape(p) ) ) stop 'r dimension error in kick'

    r = r + (dt/m) * p
    
  end subroutine drift
  
  subroutine rattle_a ( box, dt, q, r, p, iter )
    use maths_module, only : pbc
    implicit none
    real, dimension(3),   intent(in)    :: box  ! Box lengths
    real,                 intent(in)    :: dt   ! Time step
    real, dimension(:,:), intent(in)    :: q    ! Old positions (3,n)
    real, dimension(:,:), intent(inout) :: r    ! Positions (3,n)
    real, dimension(:,:), intent(inout) :: p    ! Momenta (3,n)
    integer,              intent(out)   :: iter ! Number of iterations

    ! Stage A of RATTLE, very similar to SHAKE
    ! Iteratively adjusts positions and momenta to satisfy bond constraints.

    ! On entry to this routine we assume:
    ! q stores the positions at the start of the step ("old positions")
    ! r stores the end-of-step positions following the unconstrained drift and
    ! p stores the half-step momenta following the first unconstrained half-kick
    ! All masses are the same, which slightly simplifies some of the expressions
    ! On return from this routine, r and p will hold the constrained values
    ! m and d are module parameters declared at the top of this file

    logical, dimension(size(r,dim=2)) :: move, moved ! Iterative flags (n)

    real, dimension(3) :: rij, qij, dr
    logical            :: done
    real               :: sigma, a, lambda
    integer            :: i, j, n

    real,    parameter :: tol = 1.0e-9,  tol2 = 2.0 * tol
    integer, parameter :: iter_max = 500

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
             rij   = pbc ( rij, box )   ! Periodic boundary correction
             sigma = d**2 - sum(rij**2) ! Amount by which constraint is violated

             if ( abs(sigma) > tol2*d**2 ) then ! Test whether constraint not already satisfied

                qij = q(:,i) - q(:,j)  ! Old bond vector
                qij = pbc ( qij, box ) ! Periodic boundary correction

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

  subroutine rattle_b ( box, r, p, iter )
    use maths_module, only : pbc
    implicit none
    real, dimension(3),   intent(in)    :: box  ! Box lengths
    real, dimension(:,:), intent(in)    :: r    ! Positions (3,n)
    real, dimension(:,:), intent(inout) :: p    ! Momenta (3,n)
    integer,              intent(out)   :: iter ! Number of iterations

    ! Stage B of RATTLE
    ! This subroutine iteratively adjusts the momenta
    ! to satisfy the time derivatives of the bond constraints

    ! On entry to this routine we assume:
    ! r stores the positions at the end of the step with constraints applied
    ! p stores the momenta following the second unconstrained half-kick
    ! All masses are the same, which simplifies some of the expressions
    ! On return from this routine, p will hold the constrained values
    ! m and d are module parameters declared at the top of this file

    logical, dimension(size(r,dim=2)) :: move, moved ! Iterative flags (n)

    real, dimension(3) :: rij, vij, dp
    logical            :: done
    real               :: tau, mu, b
    integer            :: i, j, n
    real,    parameter :: tol = 1.0e-9
    integer, parameter :: iter_max = 500

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
             rij = pbc ( rij, box )            ! Periodic boundary correction
             vij = p(:,i)/m - p(:,j)/m         ! Relative velocity
             tau = -m*dot_product ( rij, vij ) ! Bond derivative violation
             b   = 2*d**2                      ! B

             mu = tau / b

             if ( abs ( mu ) > tol ) then ! Test whether constraint already satisfied

                dp      = rij * mu    ! Momentum adjustment
                p(:,i)  = p(:,i) + dp ! Adjust i momentum
                p(:,j)  = p(:,j) - dp ! Adjust j momentum
                move(i) = .true.      ! Flag that we moved i
                move(j) = .true.      ! Flag that we moved j
                done    = .false.     ! Flag that we moved something

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

  subroutine milcshake_a ( box, dt, q, r, p, iter )
    use maths_module, only : pbc, tridiagonal
    implicit none
    real, dimension(3),   intent(in)    :: box  ! Box lengths
    real,                 intent(in)    :: dt   ! Time step
    real, dimension(:,:), intent(in)    :: q    ! Old positions (3,n)
    real, dimension(:,:), intent(inout) :: r    ! Positions (3,n)
    real, dimension(:,:), intent(inout) :: p    ! Momenta (3,n)
    integer,              intent(out)   :: iter ! Number of iterations

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

    integer                            :: n ! number of atoms
    integer                            :: c ! number of constraints, n-1

    real, dimension(3,size(r,dim=2)-1) :: qij ! old bond vectors (3,c)
    real, dimension(3,size(r,dim=2)-1) :: rij ! new bond vectors (3,c)
    real, dimension(3,size(r,dim=2)-1) :: dr  ! position update vectors (3,c)

    real, dimension(size(r,dim=2)-1) :: lambda ! multipliers (c)
    real, dimension(size(r,dim=2)-1) :: sigma  ! amounts by which constraints are violated (c)
    real, dimension(size(r,dim=2)-1) :: ad     ! diagonal elements of A matrix (c)
    real, dimension(size(r,dim=2)-2) :: au     ! upper-diagonal elements of A matrix (c-1)
    real, dimension(size(r,dim=2)-2) :: al     ! lower-diagonal elements of A matrix (c-1)

    real,    parameter :: tol = 1.0e-9, tol2 = 2.0 * tol
    integer, parameter :: iter_max = 500

    n = size(r,dim=2) ! Number of atoms
    c = n - 1         ! Number of constraints
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in milcshake_a'
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in milcshake_a'
    if ( any ( shape(q) /= [3,n] ) ) stop 'q dimension error in milcshake_a'

    qij = q(:,1:c) - q(:,2:n) ! Old bond vectors
    qij = pbc ( qij, box )    ! Periodic boundary correction
    rij = r(:,1:c) - r(:,2:n) ! New (unconstrained) bond vectors
    rij = pbc ( rij, box )    ! Periodic boundary correction

    ! Elements of tridiagonal matrix A (dot products of old and new bond vectors)
    al(1:c-1) = -2*sum( qij(:,1:c-1)*rij(:,2:c),   dim=1 ) ! c-1 elements of lower-diagonal
    ad(1:c)   =  4*sum( qij(:,1:c)  *rij(:,1:c),   dim=1 ) ! c elements of diagonal
    au(1:c-1) = -2*sum( qij(:,2:c)  *rij(:,1:c-1), dim=1 ) ! c-1 elements of upper-diagonal

    iter = 0

    do ! Iterative loop until done

       rij   = r(:,1:c) - r(:,2:n)      ! New bond vectors
       rij   = pbc ( rij, box )         ! Periodic boundary correction
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

  subroutine milcshake_b ( box, r, p, iter )
    use maths_module, only : pbc, tridiagonal
    implicit none
    real, dimension(3),   intent(in)    :: box  ! Box lengths
    real, dimension(:,:), intent(in)    :: r    ! Positions (3,n)
    real, dimension(:,:), intent(inout) :: p    ! Momenta (3,n)
    integer,              intent(out)   :: iter ! Number of iterations (always 1)

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

    integer                            :: n   ! number of atoms
    integer                            :: c   ! number of constraints, n-1
    real, dimension(3,size(r,dim=2)-1) :: rij ! bond vectors (3,c)
    real, dimension(3,size(r,dim=2)-1) :: vij ! new relative velocities (3,c)
    real, dimension(3,size(r,dim=2)-1) :: dp  ! momentum change vectors (3,c)

    real, dimension(size(r,dim=2)-1) :: mu  ! multipliers (c)
    real, dimension(size(r,dim=2)-1) :: tau ! amounts by which constraints are violated (c)
    real, dimension(size(r,dim=2)-1) :: bd  ! diagonal elements of B matrix (c)
    real, dimension(size(r,dim=2)-2) :: bu  ! upper-diagonal elements of B matrix (c-1)
    real, dimension(size(r,dim=2)-2) :: bl  ! lower-diagonal elements of B matrix (c-1)

    n = size(r,dim=2) ! Number of atoms
    c = n - 1         ! Number of constraints
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in milcshake_b'
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in milcshake_b'

    rij = r(:,1:c) - r(:,2:n) ! Bond vectors
    rij = pbc ( rij, box )    ! Periodic boundary correction

    ! Elements of tridiagonal matrix B (dot products of bond vectors)
    bl =  -sum ( rij(:,1:c-1)*rij(:,2:c),   dim=1 ) ! c-1 elements of lower-diagonal
    bd = 2*sum ( rij(:,1:c)  *rij(:,1:c),   dim=1 ) ! c elements of diagonal
    bu =  -sum ( rij(:,2:c)  *rij(:,1:c-1), dim=1 ) ! c-1 elements of upper-diagonal

    vij = p(:,1:c)/m - p(:,2:n)/m         ! Relative velocities
    tau = -m*sum( vij*rij, dim=1 )        ! Amounts by which constraints are violated
    mu  = tridiagonal ( bl, bd, bu, tau ) ! Solve tridiagonal system
    dp  = spread(mu,dim=1,ncopies=3)*rij  ! Resultant update vectors

    p(:,1:c) = p(:,1:c) + dp ! Adjust i momenta
    p(:,2:n) = p(:,2:n) - dp ! Adjust j momenta

    iter = 1

  end subroutine milcshake_b

  subroutine andersen ( temperature, p )
    use maths_module, only : random_normal
    implicit none
    real,                 intent(in)  :: temperature ! Temperature
    real, dimension(:,:), intent(out) :: p           ! Momenta (3,n)

    ! Implements the thermostat of HC Andersen, J Chem Phys 72, 2384 (1980)
    ! All masses are the same
    ! m is a module parameter declared at the top of this file

    ! Warning: in this simple implementation, keep fraction=1.0 if this
    ! thermostat is to be combined with constraints. For details see:
    ! EJF Peters, N Goga, HJC Berendsen, J Chem Theor Comput 10, 4208 (2014)

    integer            :: i, x, n
    real               :: std, zeta
    real, dimension(3) :: pc
    real, parameter    :: fraction = 1.0 ! By default, thermalize all the atoms

    n =  size(p,dim=2)
    if ( size(p,dim=1) /= 3 ) stop 'p dimension error in andersen'

    std = sqrt(m*temperature)

    do i = 1, n

       call random_number ( zeta )
       if ( zeta > fraction ) cycle ! Only thermalize a fraction of atoms

       do x = 1, 3
          p(x,i) = random_normal ( 0.0, std ) ! Random momentum component
       end do

    end do

    ! Subtract centre-of-mass momentum
    pc = sum ( p, dim = 2 ) / real(n)
    p  = p - spread ( pc, dim = 2, ncopies = n )

  end subroutine andersen

  function kineng ( p ) result (k )
    implicit none
    real, dimension(:,:), intent(in) :: p ! Momenta (3,n)
    real                             :: k ! Total kinetic energy

    ! Returns total kinetic energy
    ! All masses are the same
    ! m is a module parameter declared at the top of this file

    k = 0.5*sum(p**2)/m

  end function kineng

  subroutine worst_bond ( box, r, p )
    use maths_module, only : pbc
    implicit none
    real, dimension(3),   intent(in) :: box ! Box lengths
    real, dimension(:,:), intent(in) :: r   ! Positions (3,n)
    real, dimension(:,:), intent(in) :: p   ! Momenta (3,n)

    ! Writes out 2 measures relating to bond constraints between neighbouring atoms
    ! The "deviation" is |dij-d| where dij is the magnitude of the vector rij = ri-rj
    ! The "derivative" is the time derivative of this, which turns out to be |vij.rij|/dij
    ! where vij=vi-vj is the relative velocity vector.
    ! In each case, we print the worst values for the whole molecule.
    ! Both these quantities would be zero for perfectly constrained bonds,
    ! but neither of them is dimensionless, and they are just roughly indicative of what is (nearly) conserved.
    ! All masses are the same
    ! m and d are module parameters declared at the top of this file

    integer                            :: n   ! Number of atoms
    integer                            :: c   ! Number of constraints, n-1
    real, dimension(3,size(r,dim=2)-1) :: rij ! bond vectors (3,c)
    real, dimension(3,size(r,dim=2)-1) :: vij ! relative velocities (3,c)

    real, dimension(size(r,dim=2)-1) :: dij   ! bond lengths (c)
    real, dimension(size(r,dim=2)-1) :: sigma ! bond length deviations (c)
    real, dimension(size(r,dim=2)-1) :: tau   ! bond length derivatives (c)

    n = size(r,dim=2) ! Number of atoms
    c = n - 1         ! Number of constraints
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in worst_bond'
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in worst_bond'

    rij   = r(:,1:c) - r(:,2:n)           ! Relative position vectors
    rij   = pbc ( rij, box )              ! Periodic boundary correction
    dij   = norm2(rij,dim=1)              ! Distances between atoms
    sigma = abs(dij-d)                    ! Deviations from specified bond lengths
    vij   = p(:,1:c)/m - p(:,2:n)/m       ! Relative velocity vectors
    tau   = abs(sum(rij*vij,dim=1)) / dij ! Time derivatives

    write ( unit=output_unit, fmt='(a,t31,es14.4)') 'Worst bond length deviation ', maxval(sigma)
    write ( unit=output_unit, fmt='(a,t31,es14.4)') 'Worst bond length derivative', maxval(tau)

  end subroutine worst_bond

end module md_module
