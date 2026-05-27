! md_springs.f90
! Molecular dynamics, NVE or NVT ensemble, chain molecule, multiple timesteps option
program md_springs

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Takes in a configuration of atoms in a linear chain (positions, momenta)
  ! No periodic boundary conditions
  ! Conducts molecular dynamics with springs and an option to use multiple timesteps
  ! Uses no special neighbour lists

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
  use force_module,     only : force_intro, spring_intro, force, spring
  use md_module,        only : md_intro, kick, drift, worst_bond, kineng, thermostat
  use hdf5_module,      only : open_hdf5_file, close_hdf5_file, write_hdf5_dataset, write_hdf5_attribute

  implicit none

  ! Most important variables
  integer  :: n           ! Number of atoms
  integer  :: nstep       ! Number of (long) time steps
  integer  :: n_mts       ! Number of short steps per long step
  real(dp) :: dt          ! Time step (short)
  real(dp) :: u           ! Total nonbonded potential energy
  real(dp) :: v           ! Total spring harmonic potential energy
  real(dp) :: k           ! Total kinetic energy
  logical  :: nvt         ! Option for NVT ensemble
  real(dp) :: temperature ! Required temperature

  ! Coordinates
  real(dp), allocatable :: r(:,:) ! Positions (3,n)
  real(dp), allocatable :: p(:,:) ! Momenta (3,n)
  real(dp), allocatable :: f(:,:) ! Non-bonded (slow) forces (3,n)
  real(dp), allocatable :: g(:,:) ! Spring bond (fast) forces (3,n)

  ! Arrays for output datasets
  real(dp), allocatable :: u_data(:) ! Total potential energy, nonbonded (nstep)
  real(dp), allocatable :: v_data(:) ! Total potential energy, springs (nstep)
  real(dp), allocatable :: k_data(:) ! Total kinetic energy (nstep)

  integer  :: step, stp, ioerr, nfree
  real(dp) :: cpu0, cpu

  character(*), parameter :: xyz_file  = 'md_springs.xyz'  ! Configuration file (overwritten at end)
  character(*), parameter :: hdf5_file = 'md_springs.hdf5' ! HDF5 file for simulation results

  character(:), allocatable :: title

  namelist /nml/ nstep, dt, n_mts, nvt, temperature

  write ( unit=output_unit, fmt='(a)' ) 'md_springs'
  write ( unit=output_unit, fmt='(a)' ) 'Molecular dynamics, constant-NVE/NVT ensemble, chain molecule'
  write ( unit=output_unit, fmt='(a)' ) 'No periodic boundaries'
  call force_intro
  call spring_intro
  call md_intro

  call random_init ( repeatable=.false., image_distinct=.true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nstep       = 100000
  n_mts       = 1
  dt          = 0.0002_dp
  nvt         = .false.
  temperature = 1.0_dp

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  read ( unit=input_unit, nml=nml, iostat=ioerr )
  if ( ioerr /= 0 ) then
     write ( unit=error_unit, fmt='(a,i10)') 'Error reading namelist nml from standard input', ioerr
     if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
     if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
     stop 'Error in md_springs'
  end if

  ! Write out run parameters
  if ( n_mts == 1 ) then
     write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of steps', nstep
     write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Time step dt',    dt
  else
     write ( unit=output_unit, fmt='(a)'           ) 'Multiple time steps'
     write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of long steps',      nstep
     write ( unit=output_unit, fmt='(a,t31,i10  )' ) 'Short steps per long step', n_mts
     write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Long time step',            dt*n_mts
     write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Short time step',           dt
  end if
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Total simulation time',     nstep*n_mts*dt
  if ( nvt ) then
     write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'NVT ensemble, temperature', temperature
  else
     write ( unit=output_unit, fmt='(a)') 'NVE ensemble'
  end if

  ! Read in initial configuration and allocate necessary arrays

  call read_config ( xyz_file, r, p ) ! Get r and p
  n     = size(r,dim=2) ! Number of atoms
  nfree = 3*n - 6       ! Number of degrees of freedom (conserved total linear and angular momentum)
  if ( n <= 1 ) stop 'Expect n>1 in md_springs'
  allocate ( f(3,n), g(3,n) )
  write ( unit=output_unit, fmt='(a,t31,i10)'    ) 'Number of particles', n

  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Centre of mass position', sum(r,dim=2)/n
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total linear momentum',   sum(p,dim=2)
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total angular momentum',  sum(cross_products(r,p),dim=2)
  call worst_bond ( r, p )

  ! Initial calculation of nonbonded forces, spring forces, and potential energies
  call force  ( r, u, f )
  call spring ( r, v, g )

  ! Allocate arrays for output datasets
  allocate ( u_data(nstep), v_data(nstep), k_data(nstep) )

  ! Initialize counters
  call cpu_time ( cpu0 )

  ! Column headings
  write ( unit=output_unit, fmt='(2a10)' ) 'Step', 'CPU'

  do step = 1, nstep ! Begin loop over steps

     if ( nvt ) call thermostat ( temperature, nfree, p )

     if ( n_mts == 1 ) then
        ! Single time step of length dt
        call kick   ( dt/2, f+g, p ) ! Kick half-step with nonbonded and spring forces
        call drift  ( dt, p, r )     ! Drift step
        call force  ( r, u, f )      ! Evaluate nonbonded forces and potential
        call spring ( r, v, g )      ! Evaluate spring forces and potential
        call kick   ( dt/2, f+g, p ) ! Kick half-step with nonbonded and spring forces
     else
        ! Single time step of length n_mts*dt
        call kick ( n_mts*dt/2, f, p ) ! Kick half-step (long) with nonbonded forces
        do stp = 1, n_mts ! Loop over n_mts steps of length dt
           call kick   ( dt/2, g, p ) ! Kick half-step (short) with spring forces
           call drift  ( dt, p, r )   ! Drift step (short)
           call spring ( r, v, g )    ! Evaluate spring forces and potential
           call kick   ( dt/2, g, p ) ! Kick half-step (short) with spring forces
        end do ! End loop over n_mts steps of length dt
        call force ( r, u, f )          ! Evaluate nonbonded forces and potential
        call kick  ( n_mts*dt/2, f, p ) ! Kick half-step (long) with nonbonded forces
     end if

     k = kineng ( p ) ! Calculate kinetic energy

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
  call open_hdf5_file ( hdf5_file )
  title = 'Molecular dynamics, springs, '
  if ( n_mts /= 1 ) title = title // 'multiple time steps, '
  if ( nvt ) then
     title = title // 'NVT'
     call write_hdf5_attribute ( 'Title', title )
     call write_hdf5_attribute ( 'T', temperature )
  else
     title = title // 'NVE'
     call write_hdf5_attribute ( 'Title', title )
  end if
  call write_hdf5_attribute ( 'nstep', nstep )
  call write_hdf5_attribute ( 'n_mts', n_mts )
  call write_hdf5_attribute ( 'dt', dt )
  call write_hdf5_attribute ( 'N', n )
  call write_hdf5_attribute ( 'CPU', cpu-cpu0 )

  ! Write out step-by-step values as datasets
  call write_hdf5_dataset ( 'U', u_data )
  call write_hdf5_dataset ( 'V', v_data )
  call write_hdf5_dataset ( 'K', k_data )

  ! Close HDF5 file
  call close_hdf5_file

  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Centre of mass position', sum(r,dim=2)/n
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total linear momentum',   sum(p,dim=2)
  write ( unit=output_unit, fmt='(a,t31,3f10.4)' ) 'Total angular momentum',  sum(cross_products(r,p),dim=2)
  call worst_bond ( r, p )

  ! Write out the final configuration
  call write_config(xyz_file,r,p)

  deallocate ( r, p, f, g )
  deallocate ( u_data, v_data, k_data )

  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU = ', cpu-cpu0

end program md_springs

