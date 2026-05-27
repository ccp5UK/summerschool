! mc_gibbs.f90
! Monte Carlo, Gibbs ensemble
program mc_gibbs

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
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

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env,  only : input_unit, output_unit, error_unit, iostat_end, iostat_eor, dp => real64

  use maths_module,     only : trigger
  use mc_module,        only : r_move, insert, n_swap, v_swap
  use config_io_module, only : read_config, write_config
  use potential_module, only : potential_intro, pot_all
  use hdf5_module,      only : open_hdf5_file, close_hdf5_file, write_hdf5_dataset, write_hdf5_attribute

  implicit none

  ! Most important variables
  integer  :: n(2)        ! Number of atoms in each box
  real(dp) :: box(3,2)    ! Box lengths for both systems
  real(dp) :: dr_max      ! Maximum MC displacement
  real(dp) :: dv_max      ! Maximum MC volume change
  real(dp) :: temperature ! Specified temperature
  real(dp) :: beta        ! Inverse temperature

  ! Coordinates
  real(dp), allocatable :: r1(:,:) ! Temporary array of coordinates for system 1 (3,n(1))
  real(dp), allocatable :: r2(:,:) ! Temporary array of coordinates for system 2 (3,n(2))
  real(dp), allocatable :: r(:,:)  ! Combined array of coordinates (3,n(1)+n(2))

  ! Arrays for output datasets
  real(dp), allocatable :: u_data(:,:) ! Total potential energy (2,nstep)
  real(dp), allocatable :: w_data(:,:) ! Total virial (2,nstep)
  real(dp), allocatable :: z_data(:,:) ! Widom estimate (1/z) = exp(-beta*mu) (2,nw)
  real(dp), allocatable :: v_data(:,:) ! Volume (2,nstep)
  integer,  allocatable :: n_data(:,:) ! Number of atoms (2,nstep)

  integer  :: step, nstep, i, nswap, nwidom, iw, nw, n_accepted, v_accepted, ioerr
  real(dp) :: v_percent, n_percent, cpu0, cpu
  logical  :: overlap
  real(dp) :: v(2), u(2), w(2), r_percent(2)
  integer  :: r_tried(2), r_accepted(2)

  ! Old and new configuration files, and output HDF5 file
  character(*), parameter :: xyz_file(2) = ['mc_gibbs_1.xyz','mc_gibbs_2.xyz'] ! Configuration files
  character(*), parameter :: hdf5_file   = 'mc_gibbs.hdf5'                     ! HDF file for results

  namelist /nml/ nstep, nswap, nwidom, temperature, dr_max, dv_max

  write( unit=output_unit, fmt='(a)' ) 'mc_gibbs'
  write( unit=output_unit, fmt='(a)' ) 'Monte Carlo, Gibbs ensemble'
  call potential_intro

  call random_init ( repeatable=.false., image_distinct=.true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  nstep       = 20000
  nswap       = 20
  nwidom      = 20
  temperature = 0.95_dp
  dr_max      = 0.15_dp
  dv_max      = 10.0_dp

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

  ! Read in initial configurations and allocate combined position array
  call read_config ( xyz_file(1), box(:,1), r1 ) ! Get box and r for system 1
  call read_config ( xyz_file(2), box(:,2), r2 ) ! Get box and r for system 2
  n = [ size(r1,dim=2), size(r2,dim=2) ]
  if ( n(1) <= 0 .or. n(2) <= 0 ) stop 'Expect both n > 0 in mc_gibbs'
  allocate ( r(3,sum(n)) ) ! This array contains the coordinates of both systems
  r(:,:n(1)  ) = r1
  r(:,n(1)+1:) = r2
  deallocate(r1,r2)
  v = product ( box, dim=1 ) ! Volumes of both systems
  write ( unit=output_unit, fmt='(a,t31,2i10)'   ) 'Numbers of particles', n
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
  write ( unit=output_unit, fmt='(6a10)' ) 'Step', 'CPU', 'Move 1 %', 'Move 2 %', 'V swap %', 'N swap %'

  do step = 1, nstep ! Begin loop over steps

     ! Single particle moves
     r_tried = r_tried + n
     do i = 1, n(1)
        call r_move ( beta, dr_max, box(:,1), r(:,:n(1)), u(1), w(1), r_accepted(1) )
     end do
     do i = 1, n(2)
        call r_move ( beta, dr_max, box(:,2), r(:,n(1)+1:), u(2), w(2), r_accepted(2) )
     end do

     ! Volume exchange move
     call v_swap ( beta, dv_max, n, box, r, u, w, v_accepted )

     ! Particle exchange moves
     do i = 1, nswap
        call n_swap ( beta, box, n, r, u, w, n_accepted )
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
        r_percent = 100*real(r_accepted,dp)/real(r_tried,dp)    ! Moves in both systems
        v_percent = 100*real(v_accepted,dp)/real(step,dp)       ! Volume exchanges
        n_percent = 100*real(n_accepted,dp)/real(nswap*step,dp) ! Particle exchanges
        call cpu_time ( cpu )
        write ( unit=output_unit, fmt='(i10,5f10.2)' ) step, cpu-cpu0, r_percent, v_percent, n_percent
     end if

  end do ! End loop over steps

  call cpu_time ( cpu )

  ! Create new HDF5 file and write simulation parameters as attributes
  call open_hdf5_file       ( hdf5_file )
  call write_hdf5_attribute ( 'Title', 'Monte Carlo, Gibbs' )
  call write_hdf5_attribute ( 'nstep', nstep )
  call write_hdf5_attribute ( 'T', temperature )
  call write_hdf5_attribute ( 'N1+N2', sum(n) )
  call write_hdf5_attribute ( 'V1+V2', sum(v) )

  ! Write out step values as datasets
  call write_hdf5_dataset ( 'U', u_data )
  call write_hdf5_dataset ( 'W', w_data )
  call write_hdf5_dataset ( 'Z', z_data )
  call write_hdf5_dataset ( 'N', n_data )
  call write_hdf5_dataset ( 'V', v_data )

  ! Close HDF5 file
  call close_hdf5_file

  ! Write out final configurations
  call write_config ( xyz_file(1), box(:,1), r(:,:n(1)  ) )
  call write_config ( xyz_file(2), box(:,2), r(:,n(1)+1:) )

  deallocate ( r )
  deallocate ( u_data, w_data, v_data, n_data, z_data )

  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU = ', cpu-cpu0

end program mc_gibbs
