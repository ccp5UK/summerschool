! mc_muca.f90
! Monte Carlo, grand canonical and multicanonical ensemble
program mc_muca

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! ACKNOWLEDGEMENT
  ! An earlier code written by Leo Lue, on which this was based, is gratefully acknowledged.
  ! Needless to say, he is in no way responsible for any mistakes in this code.
  ! There is a huge literature on this topic, and many ways to approach the problem.
  ! Here are just a few pointers.
  ! BA Berg, J Stat Phys 82, 323 (1996)
  ! BA Berg, Nucl Phys B Proc Supp 63, 982 (1998)
  ! AD Bruce and NB Wilding, Adv Chem Phys 127, 1 (2004)
  ! NB Wilding, J Phys Cond Matt 28, 414016 (2016)
  ! K Binder et al, Amer J Phys 80, 1099 (2012)

  ! Takes in a configuration of atoms (positions), which may have zero atoms
  ! and optionally a set of weights for each value of N
  ! Conducts grand ensemble MC, or multicanonical ensemble MC,
  ! at the given temperature and volume
  ! Cuboidal periodic boundary conditions
  ! Uses no special neighbour lists

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  ! All calculations are in simulation units defined by the model
  ! For example, for Lennard-Jones, sigma = 1, epsilon = 1
  ! There is nothing here specific to Lennard-Jones
  ! The model is defined in potential_module

  ! Data output uses the HDF5 library

  ! Working precision for reals, dp, is real64 throughout
  ! Standard integers should be OK, but the value of nstep must be less than huge(nstep) and huge(h)
  use, intrinsic :: iso_fortran_env,  only : input_unit, output_unit, error_unit, iostat_end, iostat_eor, dp => real64

  use maths_module,     only : trigger, cumsum
  use mc_module,        only : r_move, n_move
  use potential_module, only : potential_intro, pot_all
  use config_io_module, only : read_config, write_config
  use hdf5_module

  implicit none

  ! Most important variables
  real(dp) :: box(3)    ! Box lengths

  ! Coordinates
  real(dp), allocatable :: r(:,:)   ! Positions (3,nmax)
  real(dp), allocatable :: tmp(:,:) ! Temporary array

  ! N-dependent variables
  integer,  allocatable :: h(:)      ! Histogram of N (0:nmax)
  real(dp), allocatable :: phi(:)    ! Weights for N-moves (0:nmax)
  real(dp), allocatable :: mu_old(:) ! Old mu(N) table (nmax)
  real(dp), allocatable :: mu_new(:) ! New mu(N) table (nmax)
  real(dp), allocatable :: mu_mix(:) ! Mixture of old and new mu(N) (nmax)
  real(dp), allocatable :: wt_old(:) ! Old statistical weights (nmax)
  real(dp), allocatable :: wt_new(:) ! New statistical weights (nmax)
  real(dp), allocatable :: wt_mix(:) ! Combined old and new statistical weights (nmax)

  character(*), parameter :: xyz_file  = 'mc_muca.xyz'  ! Configuration file (overwritten at end)
  character(*), parameter :: hdf5_file = 'mc_muca.hdf5' ! HDF file containing weights & histograms

  integer  :: nmax, n, ioerr, r_accepted, n_accepted, step, nstep
  real(dp) :: beta, mu, temperature, r_percent, n_percent
  real(dp) :: dr_max, u, w, wtsum, cpu0, cpu
  logical  :: overlap, use_weights

  namelist /nml/ temperature, mu, nmax, nstep, dr_max

  write(output_unit,'(a)') 'mc_muca'
  write(output_unit,'(a)') 'Grand canonical / multicanonical simulation of system of atoms'
  call potential_intro

  call random_init ( repeatable=.false., image_distinct=.true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  temperature = 0.95_dp
  mu          = -3.14_dp
  nmax        = 285
  nstep       = 2000000 
  dr_max      = 0.15_dp

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  read ( unit=input_unit, nml=nml, iostat=ioerr )
  if ( ioerr /= 0 ) then
     write ( unit=error_unit, fmt='(a,i10)') 'Error reading namelist nml from standard input', ioerr
     if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
     if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
     stop 'Error in muca'
  end if

  beta = 1 / temperature

  ! Write out run parameters
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of steps        ',  nstep
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Maximum number of atoms',  nmax
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Temperature T',            temperature
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Maximum displacement',     dr_max
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Chemical potential mu',    mu

  ! Allocate arrays
  allocate ( h(0:nmax), phi(0:nmax) )
  allocate ( mu_old(nmax), mu_new(nmax), mu_mix(nmax) )
  allocate ( wt_old(nmax), wt_new(nmax), wt_mix(nmax) )

  ! Read mu table and statistical weights, if the HDF5 file exists
  inquire ( file=hdf5_file, exist=use_weights )

  if ( use_weights ) then

     write ( unit=output_unit, fmt='(a)' ) 'Multicanonical simulation'

     ! Open existing HDF5 file, for reading only
     call open_hdf5_file ( hdf5_file, readonly=.true. )

     ! Read nmax from attributes and confirm that it matches, for simplicity
     call read_hdf5_attribute ( 'Nmax', n )
     if ( n /= nmax ) stop 'We insist on a matching value of Nmax in this file'

     ! Read mu table and statistical weights from datasets
     write ( unit=output_unit, fmt='(a)' ) 'Reading multicanonical mu and weights '
     call read_hdf5_dataset ( 'mu', mu_old )
     call read_hdf5_dataset ( 'wt', wt_old )

     ! Close HDF5 file
     call close_hdf5_file

  else

     ! Set mu table for standard grand canonical simulation
     ! These values will persist through successive simulations until over-written 
     write ( unit=output_unit, fmt='(a)' ) 'Grand canonical simulation'
     mu_old = mu
     wt_old = 0
  end if

  ! phi(n) = phi(n-1) + mu(n), n = 1, 2, ..., nmax
  ! Uses cumulative sum function, with phi(0)=0
  phi(0:nmax)  = [ 0.0_dp, cumsum(mu_old) ]

  ! Read initial configuration
  call read_config ( xyz_file, box, r ) ! Get box and positions
  n = size(r,dim=2)
  if ( n > nmax ) stop 'Configuration n > nmax'
  allocate ( tmp(3,nmax) )
  tmp(:,:n) = r
  call move_alloc ( tmp, r )
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Initial n', n
  write ( unit=output_unit, fmt='(a,t31,3f10.4)') 'Box',       box

  ! Check for overlaps; calculate total potential energy & virial
  call pot_all ( box, r(:,:n), u, w, overlap )
  if ( overlap ) stop 'Overlaps detected in initial configuration'

  ! Initialize histograms and counters
  call cpu_time ( cpu0 )
  h = 0
  r_accepted = 0
  n_accepted = 0

  ! Column headings
  write ( unit=output_unit, fmt='(4a10)' ) 'Step', 'CPU', 'Move %', 'N move %' 

  do step = 1, nstep ! Begin loop over steps

     call r_move ( beta, dr_max, box, r(:,:n), u, w, r_accepted )

     call n_move ( beta, phi, box, n, r, u, w, n_accepted )

     ! Accumulate histogram
     if ( n <= nmax ) h(n) = h(n) + 1

     ! Write out move acceptance rates
     if ( trigger(step) ) then
        r_percent = 100*real(r_accepted,dp)/real(step,dp)
        n_percent = 100*real(n_accepted,dp)/real(step,dp)
        call cpu_time ( cpu )
        write ( unit=output_unit, fmt='(i10,3f10.2)' ) step, cpu-cpu0, r_percent, n_percent
     end if

  end do ! End loop over steps

  call cpu_time ( cpu )

  ! Write final configuration
  call write_config ( xyz_file, box, r(:,:n) )

  ! Calculate new estimates of mu and statistical weights
  do n = 1, nmax
     mu_new(n) = 0.0_dp
     wt_new(n) = 0.0_dp
     if ( h(n-1) > 0 .and. h(n) > 0 ) then
        mu_new(n) = -temperature*log(real(h(n),dp)/real(h(n-1),dp)) + mu_old(n)
        wt_new(n) = real(h(n-1),dp)*real(h(n),dp)/(real(h(n-1),dp)+real(h(n),dp))
     end if
  end do

  ! Mix new and old mu and statistical weights
  do n = 1, nmax
     wtsum = wt_old(n) + wt_new(n)
     if ( wtsum > tiny(1.0_dp) ) then
        mu_mix(n) = (wt_old(n)*mu_old(n) + wt_new(n)*mu_new(n)) / wtsum
        wt_mix(n) = wtsum
     else
        mu_mix(n) = mu_old(n)
        wt_mix(n) = wt_old(n)
     end if
  end do

  ! Create new HDF5 file, write simulation parameters as attributes
  call open_hdf5_file ( hdf5_file )
  call write_hdf5_attribute ( 'Title', 'Monte Carlo MUCA' )
  call write_hdf5_attribute ( 'Nmax', nmax )
  call write_hdf5_attribute ( 'nstep', nstep )
  call write_hdf5_attribute ( 'L', box )
  call write_hdf5_attribute ( 'V', product(box) )
  call write_hdf5_attribute ( 'T', temperature )
  call write_hdf5_attribute ( 'mu', mu )

  ! Write histogram, mu table, and statistical weights as datasets
  call write_hdf5_dataset   ( 'h', h )
  call write_hdf5_dataset   ( 'mu', mu_mix )
  call write_hdf5_dataset   ( 'wt', wt_mix )

  ! Close HDF5 file
  call close_hdf5_file

  deallocate ( r )
  deallocate ( h, phi )
  deallocate ( mu_old, wt_old, mu_new, wt_new, mu_mix, wt_mix )
  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU = ', cpu-cpu0

end program mc_muca

