! md_constraints.f90
! Molecular dynamics, constraints, NVE or NVT ensemble, chain molecule

program md_constraints

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Takes in a configuration of atoms in a linear chain (positions, momenta)
  ! No periodic boundary conditions
  ! Conducts molecular dynamics with constraints (RATTLE or MILC SHAKE)
  ! Uses no special neighbour lists
  ! There is no periodic box

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Input configuration, all calculations, and all results 
  ! are given in simulation units defined by the model
  ! which is defined in the associated modules.
  ! For example, Lennard-Jones, sigma = 1, epsilon = 1, mass = 1

  ! Data output uses the HDF5 library

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env,  only : input_unit, output_unit, error_unit, iostat_end, iostat_eor, dp => real64

  use maths_module,     only : trigger, cross_products
  use config_io_module, only : read_config, write_config
  use force_module,     only : force_intro, force
  use md_module,        only : md_intro, constraint_intro, kick, drift, worst_bond, kineng, thermostat, &
       &                       rattle_a, rattle_b, milcshake_a, milcshake_b
  use hdf5_module,      only : open_hdf5_file, close_hdf5_file, write_hdf5_dataset, write_hdf5_attribute

  implicit none

  ! Most important variables
  integer  :: n           ! Number of atoms
  integer  :: nstep       ! Number of steps
  real(dp) :: dt          ! Time step
  real(dp) :: u           ! Potential energy
  real(dp) :: k           ! Total kinetic energy
  logical  :: nvt         ! Option for NVT ensemble
  real(dp) :: temperature ! Required temperature
  logical  :: milcshake   ! Option for constraint algorithm

  ! Coordinates
  real(dp), allocatable :: r(:,:) ! Positions (3,n)
  real(dp), allocatable :: q(:,:) ! Old positions (3,n)
  real(dp), allocatable :: p(:,:) ! Momenta (3,n)
  real(dp), allocatable :: f(:,:) ! Non-bonded forces (3,n)

  ! Arrays for output datasets
  real(dp), allocatable :: u_data(:) ! Total potential energy (nstep)
  real(dp), allocatable :: k_data(:) ! Total kinetic energy (nstep)

  ! Define procedure pointers for constraint routines
  procedure(rattle_a), pointer :: constraints_a => null()
  procedure(rattle_b), pointer :: constraints_b => null()

  integer  :: step, ioerr, iter, iter_a, iter_b, nfree
  real(dp) :: cpu0, cpu

  character(*), parameter :: xyz_file  = 'md_constraints.xyz'  ! Configuration file (overwritten at end)
  character(*), parameter :: hdf5_file = 'md_constraints.hdf5' ! HDF5 file for simulation results

  character(:), allocatable :: title

  namelist /nml/ nstep, dt, nvt, temperature, milcshake

  write ( unit=output_unit, fmt='(a)' ) 'md_constraints'
  write ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE/NVT ensemble, chain molecule'
  write ( unit=output_unit, fmt='(a)' ) 'No periodic boundaries'
  write ( unit=output_unit, fmt='(a)' ) 'Bond constraints'
  call force_intro
  call md_intro
  call constraint_intro

  call random_init ( repeatable=.false., image_distinct=.true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nstep       = 10000
  dt          = 0.002_dp
  nvt         = .false.
  temperature = 1.0_dp
  milcshake   = .false.

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  read ( unit=input_unit, nml=nml, iostat=ioerr )
  if ( ioerr /= 0 ) then
     write ( unit=error_unit, fmt='(a,i10)') 'Error reading namelist nml from standard input', ioerr
     if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
     if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
     stop 'Error in md_constraints'
  end if

  ! Write out run parameters
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of steps',       nstep
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Time step dt',          dt
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Total simulation time', nstep*dt
  if ( nvt ) then
     write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'NVT ensemble, temperature', temperature
  else
     write ( unit=output_unit, fmt='(a)'   ) 'NVE ensemble'
  end if

  ! Read in initial configuration and allocate necessary arrays

  call read_config ( xyz_file, r, p ) ! Get r and p
  n     = size(r,dim=2)   ! Number of atoms
  nfree = 3*n - 6 - (n-1) ! Number of degrees of freedom (conserved total linear and angular momentum, and bonds)
  if ( n <= 1 ) stop 'Expect n>1 in md_constraints'
  allocate ( q(3,n), f(3,n) )
  write ( unit=output_unit, fmt='(a,t31,i10)'    ) 'Number of particles', n
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Centre of mass position', sum(r,dim=2)/n
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total linear momentum',   sum(p,dim=2)
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total angular momentum',  sum(cross_products(r,p),dim=2)
  call worst_bond ( r, p )

  ! Calculate initial values and forces
  call force ( r, u, f )

  ! Allocate arrays for output datasets
  allocate ( u_data(nstep), k_data(nstep) )

  ! Point to selected constraint routines
  if ( milcshake ) then
     write ( unit=output_unit, fmt='(a)' ) 'MILC SHAKE constraint method'
     constraints_a => milcshake_a
     constraints_b => milcshake_b
  else
     write ( unit=output_unit, fmt='(a)' ) 'RATTLE constraint method'
     constraints_a => rattle_a
     constraints_b => rattle_b
  end if

  ! Initialize counters
  call cpu_time ( cpu0 )
  iter_a = 0
  iter_b = 0

  ! Column headings
  write ( unit=output_unit, fmt='(4a10)' ) 'Step', 'CPU', 'Iter A', 'Iter B'

  do step = 1, nstep ! Begin loop over steps

     if ( nvt ) call thermostat ( temperature, nfree, p )

     q = r                                    ! Old positions needed in constraints_a
     call kick ( dt/2, f, p )                 ! Kick half-step
     call drift ( dt, p, r )                  ! Drift step
     call constraints_a ( dt, q, r, p, iter ) ! Bond length constraints
     iter_a = iter_a + iter                   ! Accumulate A iterations
     call force ( r, u, f )                   ! Force evaluation
     call kick ( dt/2, f, p )                 ! Kick half-step
     call constraints_b ( r, p, iter )        ! Bond length derivative constraints
     iter_b = iter_b + iter                   ! Accumulate B iterations

     k = kineng ( p ) ! Calculate kinetic energy

     ! Save data for this step
     u_data(step) = u ! Potential energy
     k_data(step) = k ! Kinetic energy

     ! Write out step counter and average iteration counts
     if ( trigger(step) ) then
        call cpu_time ( cpu )
        write ( unit=output_unit, fmt='(i10,3f10.2)' ) step, cpu-cpu0, real(iter_a,dp)/step, real(iter_b,dp)/step
     end if

  end do ! End loop over steps

  call cpu_time ( cpu )

  ! Create new HDF5 file and write simulation parameters as attributes
  call open_hdf5_file ( hdf5_file )
  title = 'Molecular dynamics, constraints, '
  if ( nvt ) then
     title = title // 'NVT'
     call write_hdf5_attribute ( 'Title', title )
     call write_hdf5_attribute ( 'T', temperature )
  else
     title = title // 'NVE'
     call write_hdf5_attribute ( 'Title', title )
  end if
  call write_hdf5_attribute ( 'nstep', nstep )
  call write_hdf5_attribute ( 'dt', dt )
  call write_hdf5_attribute ( 'N', n )
  call write_hdf5_attribute ( 'CPU', cpu-cpu0 )

  ! Write out step-by-step values as datasets
  call write_hdf5_dataset ( 'U', u_data )
  call write_hdf5_dataset ( 'K', k_data )

  ! Close HDF5 file
  call close_hdf5_file

  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Centre of mass position', sum(r,dim=2)/n
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total linear momentum',   sum(p,dim=2)
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total angular momentum',  sum(cross_products(r,p),dim=2)
  call worst_bond ( r, p )

  ! Write out the final configuration
  call write_config(xyz_file,r,p)

  deallocate ( r, q, p, f )
  deallocate ( u_data, k_data )

  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU = ', cpu-cpu0

end program md_constraints

