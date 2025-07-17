! mc_gibbs.f90
! Monte Carlo, Gibbs ensemble
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
program mc_gibbs

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Takes in a pair of configurations of atoms (positions)
  ! Cuboidal periodic boundary conditions
  ! Conducts Gibbs ensemble Monte Carlo at the given temperature T,
  ! keeping total volume V1+V2 and total number of atoms N1+N2 fixed.
  ! To avoid some inconvenient tests, we disallow configurations in which either box is empty
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! All calculations are in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1
  ! There is nothing here specific to Lennard-Jones
  ! The model is defined in potential_module

  ! Data output uses the HDF5 library

  use, intrinsic :: iso_fortran_env,  only : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  use               maths_module,     only : trigger
  use               mc_module,        only : r_move, insert
  use               mc_gibbs_module,  only : n_move, v_move
  use               config_io_module, only : n_config, read_config, write_config
  use               potential_module, only : introduction, pot_all
  use               hdf5_module,      only : open_file, close_file, write_dataset, write_attribute

  implicit none

  ! Most important variables
  integer, dimension(2)                :: n           ! Number of atoms in each box
  real,    dimension(3,2)              :: box         ! Box lengths for both systems
  real,    dimension(:,:), allocatable :: r           ! Positions for both systems (3,n(1)+n(2))
  real                                 :: dr_max      ! Maximum MC displacement
  real                                 :: dv_max      ! Maximum MC volume change
  real                                 :: temperature ! Specified temperature
  real                                 :: beta        ! Inverse temperature

  ! Arrays for output datasets
  real,    dimension(:,:), allocatable :: u_data ! Total potential energy (2,nstep)
  real,    dimension(:,:), allocatable :: w_data ! Total virial (2,nstep)
  real,    dimension(:,:), allocatable :: z_data ! Widom estimate (1/z) = exp(-beta*mu) (2,nw)
  real,    dimension(:,:), allocatable :: v_data ! Volume (2,nstep)
  integer, dimension(:,:), allocatable :: n_data ! Number of atoms (2,nstep)

  integer               :: step, nstep, i, nswap, nwidom, iw, nw, n_accepted, v_accepted, ioerr
  real                  :: v_percent, n_percent, cpu0, cpu
  logical               :: overlap
  real,    dimension(2) :: v, u, w, r_percent
  integer, dimension(2) :: r_tried, r_accepted

  ! Old and new configuration files, and output HDF5 file
  character(len=*), dimension(2), parameter :: config_old = ['config_one_old.dat','config_two_old.dat']
  character(len=*), dimension(2), parameter :: config_new = ['config_one.dat','config_two.dat']
  character(len=*),               parameter :: hdf5_file   = 'mc_gibbs.hdf5'

  namelist /nml/ nstep, nswap, nwidom, temperature, dr_max, dv_max

  write( unit=output_unit, fmt='(a)' ) 'mc_gibbs'
  write( unit=output_unit, fmt='(a)' ) 'Monte Carlo, Gibbs ensemble'
  call introduction

  call random_init ( .false., .true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nstep       = 20000
  nswap       = 20
  nwidom      = 20
  temperature = 0.95
  dr_max      = 0.15
  dv_max      = 10.0

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  read ( unit=input_unit, nml=nml, iostat=ioerr )
  if ( ioerr /= 0 ) then
     write ( unit=error_unit, fmt='(a,i10)') 'Error reading namelist nml from standard input', ioerr
     if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
     if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
     stop 'Error in mc_gibbs'
  end if

  beta = 1 / temperature

  ! Write out run parameters
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of steps',         nstep
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Widom attempts per step', nwidom
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Swap attempts per step',  nswap
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Temperature',             temperature
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Maximum displacement',    dr_max
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Maximum volume change',   dv_max

  ! Read in initial configurations and allocate position array

  n(1) = n_config ( config_old(1) ) ! Get N for system 1
  n(2) = n_config ( config_old(2) ) ! Get N for system 2
  write ( unit=output_unit, fmt='(a,t31,2i10)' ) 'Numbers of particles', n
  if ( sum(n) <= 0 ) stop 'Expect total n > 0 in mc_gibbs'
  allocate ( r(3,sum(n)) ) ! This array contains the coordinates of both systems
  call read_config ( config_old(1), box(:,1), r(:,:n(1)  ) ) ! Get box and r for system 1
  call read_config ( config_old(2), box(:,2), r(:,n(1)+1:) ) ! Get box and r for system 2
  v = product ( box, dim=1 ) ! Volumes of both systems
  write ( unit=output_unit, fmt='(a,t31,2f10.4)' ) 'Box volumes', v
  write ( unit=output_unit, fmt='(a,t31,2f10.4)' ) 'Densities',   n / v

  ! Initial energy and overlap check
  call pot_all ( box(:,1), r(:,:n(1)), u(1), w(1), overlap ) 
  if ( overlap ) then
     write ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration 1'
     stop 'Error in mc_gibbs'
  end if
  call pot_all ( box(:,2), r(:,n(1)+1:), u(2), w(2), overlap ) 
  if ( overlap ) then
     write ( unit=error_unit, fmt='(a)') 'Overlap in initial configuration 2'
     stop 'Error in mc_gibbs'
  end if

  ! Allocate arrays for output datasets
  allocate ( u_data(2,nstep), w_data(2,nstep), v_data(2,nstep), n_data(2,nstep) )
  nw = nstep * nwidom
  allocate ( z_data(2,nw) )

  ! Initialize counters
  call cpu_time ( cpu0 )
  r_tried    = 0
  r_accepted = 0
  v_accepted = 0
  n_accepted = 0
  iw = 0

  ! Column headings
  write ( unit=output_unit, fmt='(6a10)' ) 'Step', 'CPU', 'Move 1 %', 'Move 2 %', 'V move %', 'N move %'

  do step = 1, nstep ! Begin loop over steps

     ! Single particle moves
     r_tried = r_tried + n
     do i = 1, n(1)
        call r_move ( beta, dr_max, box(:,1), r(:,:n(1)), u(1), w(1), r_accepted(1) )
     end do
     do i = 1, n(2)
        call r_move ( beta, dr_max, box(:,2), r(:,n(1)+1:), u(2), w(2), r_accepted(2) )
     end do

     ! Volume move
     call v_move ( beta, dv_max, n, box, r, u, w, v_accepted )

     ! Particle exchange moves
     do i = 1, nswap
        call n_move ( beta, box, n, r, u, w, n_accepted )
     end do

     ! Store step values for each box
     u_data(:,step) = u
     w_data(:,step) = w
     v_data(:,step) = product ( box, dim=1 )
     n_data(:,step) = n

     ! Widom test particle insertions
     do i = 1, nwidom
        iw = iw + 1
        if ( iw > nw ) stop 'This should be impossible'
        z_data(1,iw) = insert ( beta, box(:,1), r(:,:n(1)  ) )
        z_data(2,iw) = insert ( beta, box(:,2), r(:,n(1)+1:) )
     end do

     ! Write out move acceptance rates
     if ( trigger(step) ) then
        r_percent = 100*real(r_accepted)/real(r_tried)    ! Moves in both systems
        v_percent = 100*real(v_accepted)/real(step)       ! Volume exchanges
        n_percent = 100*real(n_accepted)/real(nswap*step) ! Particle exchanges
        call cpu_time ( cpu )
        write ( unit=output_unit, fmt='(i10,5f10.2)' ) step, cpu-cpu0, r_percent, v_percent, n_percent
     end if

  end do ! End loop over steps

  call cpu_time ( cpu )

  ! Create new HDF5 file and write simulation parameters as attributes
  call open_file       ( hdf5_file )
  call write_attribute ( 'Title', 'Monte Carlo, Gibbs' )
  call write_attribute ( 'nstep', nstep )
  call write_attribute ( 'T', temperature )
  call write_attribute ( 'N1+N2', sum(n) )
  call write_attribute ( 'V1+V2', sum(v) )

  ! Write out step values as datasets
  call write_dataset ( 'U', u_data )
  call write_dataset ( 'W', w_data )
  call write_dataset ( 'Z', z_data )
  call write_dataset ( 'N', n_data )
  call write_dataset ( 'V', v_data )

  ! Close HDF5 file
  call close_file

  ! Write out final configurations
  call write_config ( config_new(1), box(:,1), r(:,:n(1)  ) )
  call write_config ( config_new(2), box(:,2), r(:,n(1)+1:) )

  deallocate ( r )
  deallocate ( u_data, w_data, v_data, n_data, z_data )

  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU = ', cpu-cpu0

end program mc_gibbs
