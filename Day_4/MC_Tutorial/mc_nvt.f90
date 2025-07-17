! mc_nvt.f90
! Monte Carlo, constant-NVT ensemble
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
program mc_nvt

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Takes in a configuration of atoms (positions)
  ! Cuboidal periodic boundary conditions
  ! Conducts Monte Carlo at the given temperature T and volume V,
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! All calculations are in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1
  ! There is nothing here specific to Lennard-Jones
  ! The model is defined in potential_module

  ! Each step includes N single-atom moves, choosing the atoms randomly.
  ! Choosing the atoms sequentially is also correct,
  ! even though it violates the strict detailed balance condition.
  ! See for example VI Manousiouthakis and MW Deem, J Chem Phys 110, 2753 (1999).

  ! Data output uses the HDF5 library

  use, intrinsic :: iso_fortran_env,  only : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  use               maths_module,     only : trigger
  use               mc_module,        only : r_move, insert
  use               config_io_module, only : n_config, read_config, write_config
  use               potential_module, only : introduction, pot_all, t_config
  use               hdf5_module,      only : open_file, close_file, write_dataset, write_attribute

  implicit none

  ! Most important variables
  integer                             :: n           ! Number of atoms
  real,   dimension(3)                :: box         ! Box lengths
  real,   dimension(:,:), allocatable :: r           ! Positions (3,n)
  real                                :: dr_max      ! Maximum MC displacement
  real                                :: temperature ! Specified temperature
  real                                :: beta        ! Inverse temperature

  ! Arrays for output datasets
  real, dimension(:),     allocatable :: u_data ! Total potential energy (nstep)
  real, dimension(:),     allocatable :: w_data ! Total virial (nstep)
  real, dimension(:),     allocatable :: t_data ! Configurational temperature (nstep)
  real, dimension(:),     allocatable :: z_data ! Widom estimate of (1/z) = exp(-beta*mu) (nw)
  real, dimension(:,:,:), allocatable :: r_data ! Configurations at each step (3,n,nr)

  integer :: i, step, nstep, gap, nwidom, r_accepted, ir, nr, iw, nw, ioerr
  real    :: vol, u, w, r_percent, cpu0, cpu
  logical :: overlap

  character(len=*), parameter :: config_old = 'config_old.dat' ! Input configuration file
  character(len=*), parameter :: config_new = 'config.dat'     ! Output configuration file
  character(len=*), parameter :: hdf5_file  = 'mc_nvt.hdf5'    ! HDF5 file for simulation results

  namelist /nml/ nstep, gap, nwidom, temperature, dr_max

  write( unit=output_unit, fmt='(a)' ) 'mc_nvt'
  write( unit=output_unit, fmt='(a)' ) 'Monte Carlo, constant-NVT ensemble'
  call introduction

  call random_init ( .false., .true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nstep       = 20000
  gap         = 20
  nwidom      = 20
  temperature = 2.0
  dr_max      = 0.2

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  read ( unit=input_unit, nml=nml, iostat=ioerr )
  if ( ioerr /= 0 ) then
     write ( unit=error_unit, fmt='(a,i10)') 'Error reading namelist nml from standard input', ioerr
     if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
     if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
     stop 'Error in mc_nvt'
  end if

  beta = 1 / temperature

  ! Write out run parameters
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of steps',             nstep
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Gap between r saves',         gap
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Insertion attempts per step', nwidom
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Temperature',                 temperature
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Maximum displacement',        dr_max
  
  ! Read in initial configuration and allocate position array

  n = n_config ( config_old ) ! Get n
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of particles', n
  if ( n <= 0 ) stop 'Expect n>0 in mc_nvt'
  allocate ( r(3,n) )
  call read_config ( config_old, box, r ) ! Get box and r
  vol = product ( box )
  write ( unit=output_unit, fmt='(a,t31,3f10.4)') 'Box',        box
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Box volume', vol
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Density',    n / vol

  ! Initial energy and overlap check
  call pot_all ( box, r, u, w, overlap ) 
  if ( overlap ) then
     write ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration'
     stop 'Error in mc_nvt'
  end if

  ! Allocate arrays for output datasets
  allocate ( u_data(nstep), w_data(nstep), t_data(nstep) )
  nw = nstep * nwidom ! Number of test-particle insertions
  allocate ( z_data(nw) )
  nr = nstep / gap  ! Number of configurations to be stored
  allocate ( r_data(3,n,nr) )

  ! Initialize counters
  call cpu_time ( cpu0 )
  r_accepted = 0
  ir = 0
  iw = 0

  ! Column headings
  write ( unit=output_unit, fmt='(3a10)' ) 'Step', 'CPU', 'Move %'
  
  do step = 1, nstep ! Begin loop over steps

     ! N single particle moves
     do i = 1, n
        call r_move ( beta, dr_max, box, r, u, w, r_accepted )
     end do

     ! Save data for this step
     u_data(step) = u                   ! Potential energy
     w_data(step) = w                   ! Total virial
     t_data(step) = t_config ( box, r ) ! Configurational temperature

     ! Widom test particle insertions
     do i = 1, nwidom
        iw = iw + 1
        if ( iw > nw ) stop 'z_data overflow' ! Should never happen
        z_data(iw) = insert ( beta, box, r )
     end do

     ! Configuration 
     if ( mod(step,gap) == 0 ) then
        ir = ir + 1
        if ( ir > nr ) stop 'r_data overflow' ! Should never happen
        r_data(:,:,ir) = r
     end if

     ! Write out step number and move acceptance rate
     if ( trigger(step) ) then
        r_percent = 100*real(r_accepted)/real(n*step)
        call cpu_time ( cpu )
        write ( unit=output_unit, fmt='(i10,2f10.2)' ) step, cpu-cpu0, r_percent
     end if
     
  end do ! End loop over steps

  call cpu_time ( cpu )

  ! Create new HDF5 file and write simulation parameters as attributes
  call open_file       ( hdf5_file )
  call write_attribute ( 'Title', 'Monte Carlo, NVT' )
  call write_attribute ( 'nstep', nstep )
  call write_attribute ( 'T', temperature )
  call write_attribute ( 'L', box )
  call write_attribute ( 'V', vol )
  call write_attribute ( 'N', n )

  ! Write out step values as dataset
  call write_dataset ( 'U', u_data )
  call write_dataset ( 'W', w_data )
  call write_dataset ( 'Z', z_data )
  call write_dataset ( 'T', t_data )
  call write_dataset ( 'r', r_data   )

  ! Close HDF5 file
  call close_file
  
  ! Write out final configuration
  call write_config ( config_new, box, r )

  deallocate ( r )
  deallocate ( u_data, w_data, z_data, t_data, r_data )
  
  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU = ', cpu-cpu0

end program mc_nvt
