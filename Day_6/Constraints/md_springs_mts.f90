! md_springs_mts.f90
! Molecular dynamics, multiple timesteps, NVE or NVT ensemble, chain molecule
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
program md_springs_mts

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Takes in a configuration of atoms in a linear chain (positions, momenta)
  ! Periodic boundary conditions
  ! Conducts molecular dynamics with springs and multiple timesteps
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
  use               force_module,     only : force_intro, force, spring
  use               md_module,        only : md_intro, kick, drift, worst_bond, kineng, andersen
  use               hdf5_module,      only : open_file, close_file, write_dataset, write_attribute

  implicit none

  ! Most important variables
  integer                              :: n   ! Number of atoms
  real,    dimension(3)                :: box ! Box lengths
  real,    dimension(:,:), allocatable :: r   ! Positions (3,n)
  real,    dimension(:,:), allocatable :: p   ! Momenta (3,n)
  real,    dimension(:,:), allocatable :: f   ! Non-bonded (slow) forces (3,n)
  real,    dimension(:,:), allocatable :: g   ! Spring bond (fast) forces (3,n)

  ! Arrays for output datasets
  real, dimension(:), allocatable :: u_data ! Total potential energy, nonbonded (nstep)
  real, dimension(:), allocatable :: v_data ! Total potential energy, springs (nstep)
  real, dimension(:), allocatable :: k_data ! Total kinetic energy (nstep)

  integer :: nstep       ! Number of long time steps
  integer :: n_mts       ! Number of short steps per long step
  real    :: dt          ! Time step (short)
  real    :: u           ! Total nonbonded potential energy
  real    :: v           ! Total spring harmonic potential energy
  real    :: k           ! Total kinetic energy
  logical :: nvt         ! Option for NVT ensemble
  real    :: temperature ! Required temperature

  integer :: step, stp, ioerr
  real    :: cpu0, cpu

  character(len=*), parameter :: conf_file = 'md_springs.dat'      ! Input configuration file
  character(len=*), parameter :: hdf5_file = 'md_springs_mts.hdf5' ! HDF5 file for simulation results

  namelist /nml/ nstep, dt, n_mts, nvt, temperature

  write ( unit=output_unit, fmt='(a)' ) 'md_springs_mts'
  write ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE/NVT ensemble, chain molecule'
  write ( unit=output_unit, fmt='(a)' ) 'Periodic boundaries'
  write ( unit=output_unit, fmt='(a)' ) 'Multiple time steps'
  call force_intro ( springs = .true. )
  call md_intro    ( springs = .true. )

  call random_init ( .false., .true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nstep       = 2000
  n_mts       = 10
  dt          = 0.0002
  nvt         = .false.
  temperature = 1.0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  read ( unit=input_unit, nml=nml, iostat=ioerr )
  if ( ioerr /= 0 ) then
     write ( unit=error_unit, fmt='(a,i10)') 'Error reading namelist nml from standard input', ioerr
     if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
     if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
     stop 'Error in md_springs_mts'
  end if

  ! Write out run parameters
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of long steps',      nstep
  write ( unit=output_unit, fmt='(a,t31,i10  )' ) 'Short steps per long step', n_mts
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Long time step',            dt*n_mts
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Short time step',           dt
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Total simulation time',     nstep*n_mts*dt
  if ( nvt ) then
     write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'NVT ensemble, temperature', temperature
  else
     write ( unit=output_unit, fmt='(a)') 'NVE ensemble'
  end if

  ! Read in initial configuration and allocate necessary arrays

  n = n_config ( conf_file ) ! Get n
  write ( unit=output_unit, fmt='(a,t31,i10)'    ) 'Number of particles', n
  if ( n <= 0 ) stop 'Expect n>0 in md_springs_mts'
  allocate ( r(3,n), p(3,n), f(3,n), g(3,n) )
  call read_config ( conf_file, box, r, p ) ! Get r and p
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Box lengths',    box
  write ( unit=output_unit, fmt='(a,t31,f10.4)'  ) 'Density',        real(n)/product(box)
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total momentum', sum(p,dim=2)
  call worst_bond ( box, r, p )

  ! Initial calculation of nonbonded forces, spring forces, and potential energies
  call force  ( box, r, u, f )
  call spring ( box, r, v, g )

  ! Allocate arrays for output datasets
  allocate ( u_data(nstep), v_data(nstep), k_data(nstep) )

  ! Initialize counters
  call cpu_time ( cpu0 )

  ! Column headings
  write ( unit=output_unit, fmt='(2a10)' ) 'Step', 'CPU/sec'

  do step = 1, nstep ! Begin loop over steps

     if ( nvt ) call andersen ( temperature, p )

     ! Single time step of length n_mts*dt

     call kick ( 0.5*n_mts*dt, f, p ) ! Kick half-step (long) with nonbonded forces

     do stp = 1, n_mts ! Loop over n_mts steps of length dt
        call kick ( 0.5*dt, g, p )   ! Kick half-step (short) with spring forces
        call drift ( dt, p, r )      ! Drift step (short)
        call spring ( box, r, v, g ) ! Evaluate spring forces and potential
        call kick ( 0.5*dt, g, p )   ! Kick half-step (short) with spring forces
     end do ! End loop over n_mts steps of length dt

     call force ( box, r, u, f )      ! Evaluate nonbonded forces and potential
     call kick ( 0.5*n_mts*dt, f, p ) ! Kick half-step (long) with nonbonded forces

     k = kineng ( p ) ! Calculate kinetic energy
     
     ! End single time step of length n_mts*dt

     ! Save data for this step
     u_data(step) = u ! Potential energy (nonbonded part)
     v_data(step) = v ! Potential energy (spring bond part)
     k_data(step) = k ! Kinetic energy

     ! Write out step counter
     if ( trigger(step) ) then
        call cpu_time ( cpu )
        write ( unit=output_unit, fmt='(i10,f10.2)' ) step, cpu-cpu0
     end if

  end do ! End loop over steps

  call cpu_time ( cpu )

  ! Create new HDF5 file and write simulation parameters as attributes
  call open_file ( hdf5_file )
  if ( nvt ) then
     call write_attribute ( 'Title', 'Molecular dynamics multiple time steps, NVT' )
     call write_attribute ( 'T', temperature )
  else
     call write_attribute ( 'Title', 'Molecular dynamics multiple time steps, NVE' )
  end if
  call write_attribute ( 'nstep', nstep )
  call write_attribute ( 'n_mts', n_mts )
  call write_attribute ( 'dt', dt )
  call write_attribute ( 'L', box )
  call write_attribute ( 'N', n )
  call write_attribute ( 'CPU', cpu-cpu0 )

  ! Write out step-by-step values as datasets
  call write_dataset ( 'U', u_data )
  call write_dataset ( 'V', v_data )
  call write_dataset ( 'K', k_data )

  ! Close HDF5 file
  call close_file

  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total momentum', sum(p,dim=2)
  call worst_bond ( box, r, p )

  ! In this example, there is no need to write out the final configuration

  deallocate ( r, p, f, g )
  deallocate ( u_data, v_data, k_data )
  
  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU/sec = ', cpu-cpu0

end program md_springs_mts

