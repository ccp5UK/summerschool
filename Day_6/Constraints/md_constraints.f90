! md_constraints.f90
! Molecular dynamics, constraints, NVE or NVT ensemble, chain molecule
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables

! In the workshop you were asked to insert a call to constraints_b
! after the Andersen thermostat in the main loop.
! In this file the thermostatting is INCORRECT, and you need to fix it

program md_constraints

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Takes in a configuration of atoms in a linear chain (positions, momenta)
  ! Periodic boundary conditions
  ! Conducts molecular dynamics with constraints (RATTLE or MILC SHAKE)
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! Input configuration, all calculations, and all results 
  ! are given in simulation units defined by the model
  ! which is defined in the associated modules.
  ! For example, Lennard-Jones WCA, sigma = 1, epsilon = 1, mass = 1

  ! Data output uses the HDF5 library

  use, intrinsic :: iso_fortran_env,  only : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  use               maths_module,     only : trigger
  use               config_io_module, only : n_config, read_config
  use               force_module,     only : force_intro, force
  use               md_module,        only : md_intro, kick, drift, worst_bond, kineng, andersen, &
       &                                     rattle_a, rattle_b, milcshake_a, milcshake_b
  use               hdf5_module,      only : open_file, close_file, write_dataset, write_attribute

  implicit none

  ! Most important variables
  integer                             :: n     ! Number of atoms
  real,   dimension(3)                :: box   ! Box lengths
  real,   dimension(:,:), allocatable :: r     ! Positions (3,n)
  real,   dimension(:,:), allocatable :: q     ! Old positions (3,n)
  real,   dimension(:,:), allocatable :: p     ! Momenta (3,n)
  real,   dimension(:,:), allocatable :: f     ! Non-bonded forces (3,n)

  ! Arrays for output datasets
  real, dimension(:), allocatable :: u_data ! Total potential energy (nstep)
  real, dimension(:), allocatable :: k_data ! Total kinetic energy (nstep)

  ! Define procedure pointers for constraint routines
  procedure(rattle_a), pointer :: constraints_a => null()
  procedure(rattle_b), pointer :: constraints_b => null()

  integer :: nstep       ! Number of steps
  real    :: dt          ! Time step
  real    :: u           ! Potential energy
  real    :: k           ! Total kinetic energy
  logical :: nvt         ! Option for NVT ensemble
  real    :: temperature ! Required temperature
  logical :: milcshake   ! Option for constraint algorithm

  integer :: step, ioerr, iter, iter_a, iter_b
  real    :: cpu0, cpu

  character(len=*), parameter :: conf_file = 'md_constraints.dat'  ! Input configuration file
  character(len=*), parameter :: hdf5_file = 'md_constraints.hdf5' ! HDF5 file for simulation results

  namelist /nml/ nstep, dt, nvt, temperature, milcshake

  write ( unit=output_unit, fmt='(a)' ) 'md_constraints'
  write ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE/NVT ensemble, chain molecule'
  write ( unit=output_unit, fmt='(a)' ) 'Periodic boundaries'
  write ( unit=output_unit, fmt='(a)' ) 'Bond constraints'
  call force_intro ( springs = .false. )
  call md_intro    ( springs = .false. )

  call random_init ( .false., .true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nstep       = 2000
  dt          = 0.002
  nvt         = .false.
  temperature = 1.0
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

  n = n_config ( conf_file ) ! Get n
  write ( unit=output_unit, fmt='(a,t31,i10)'    ) 'Number of particles', n
  if ( n <= 0 ) stop 'Expect n>0 in md_constraints'
  allocate ( r(3,n), q(3,n), p(3,n), f(3,n) )
  call read_config ( conf_file, box, r, p ) ! Get box, r and p
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Box lengths',    box
  write ( unit=output_unit, fmt='(a,t31,f10.4)'  ) 'Density',        real(n)/product(box)
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total momentum', sum(p,dim=2)
  call worst_bond ( box, r, p )

  ! Calculate initial values and forces
  call force ( box, r, u, f )

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

     ! This thermostat will violate the constraints without the appropriate correction
     if ( nvt ) then
        call andersen ( temperature, p )
     end if

     q = r                                         ! Old positions needed in constraints_a
     call kick ( 0.5*dt, f, p )                    ! Kick half-step
     call drift ( dt, p, r )                       ! Drift step
     call constraints_a ( box, dt, q, r, p, iter ) ! Bond constraints part A
     iter_a = iter_a + iter                        ! Accumulate A iterations
     call force ( box, r, u, f )                   ! Force evaluation
     call kick ( 0.5*dt, f, p )                    ! Kick half-step
     call constraints_b ( box, r, p, iter )        ! Bond constraints part B
     iter_b = iter_b + iter                        ! Accumulate B iterations

     k = kineng ( p ) ! Calculate kinetic energy

     ! Save data for this step
     u_data(step) = u ! Potential energy
     k_data(step) = k ! Kinetic energy

     ! Write out step counter and average iteration counts
     if ( trigger(step) ) then
        call cpu_time ( cpu )
        write ( unit=output_unit, fmt='(i10,3f10.2)' ) step, cpu-cpu0, real(iter_a)/step, real(iter_b)/step
     end if

  end do ! End loop over steps

  call cpu_time ( cpu )

  ! Create new HDF5 file and write simulation parameters as attributes
  call open_file ( hdf5_file )
  if ( nvt ) then
     call write_attribute ( 'Title', 'Molecular dynamics with constraints, NVT' )
     call write_attribute ( 'T', temperature )
  else
     call write_attribute ( 'Title', 'Molecular dynamics with constraints, NVE' )
  end if
  call write_attribute ( 'nstep', nstep )
  call write_attribute ( 'dt', dt )
  call write_attribute ( 'L', box )
  call write_attribute ( 'N', n )
  call write_attribute ( 'CPU', cpu-cpu0 )

  ! Write out step-by-step values as datasets
  call write_dataset ( 'U', u_data )
  call write_dataset ( 'K', k_data )

  ! Close HDF5 file
  call close_file

  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total momentum', sum(p,dim=2)
  call worst_bond ( box, r, p )

  ! In this example, there is no need to write out the final configuration

  deallocate ( r, q, p, f )
  deallocate ( u_data, k_data )

  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU = ', cpu-cpu0

end program md_constraints

