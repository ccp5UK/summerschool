! mc_muca.f90
! Monte Carlo, grand canonical and multicanonical ensemble
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
program mc_muca

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
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

  use, intrinsic :: iso_fortran_env,  only : input_unit, output_unit, error_unit, iostat_end, iostat_eor
  use               maths_module,     only : trigger
  use               mc_module,        only : r_move
  use               mc_muca_module,   only : n_move
  use               potential_module, only : pot_all, introduction
  use               config_io_module, only : n_config, read_config, write_config
  use               hdf5_module,      only : open_file, close_file, &
       &                                     read_dataset, write_dataset, &
       &                                     read_attribute, write_attribute

  implicit none

  ! Most important variables
  real, dimension(3)                :: box    ! Box lengths
  real, dimension(:,:), allocatable :: r      ! Positions (3,nmax)
  real, dimension(:),   allocatable :: h      ! Histogram of N (0:nmax)
  real, dimension(:),   allocatable :: phi    ! Weights for N-moves (0:nmax)
  real, dimension(:),   allocatable :: mu_old ! Old mu(N) table (nmax)
  real, dimension(:),   allocatable :: mu_new ! New mu(N) table (nmax)
  real, dimension(:),   allocatable :: wt_old ! Old statistical weights (nmax)
  real, dimension(:),   allocatable :: wt_new ! New statistical weights (nmax)

  ! Old and new configuration files, old and new HDF5 files
  character(len=*), parameter :: config_old = 'config_old.dat'
  character(len=*), parameter :: config_new = 'config.dat'
  character(len=*), parameter :: hdf5_old   = 'mc_muca_old.hdf5'
  character(len=*), parameter :: hdf5_new   = 'mc_muca.hdf5'

  integer :: nmax, n, ioerr, r_accepted, n_accepted, step, nstep
  real    :: beta, mu, temperature, r_percent, n_percent
  real    :: dr_max, u, w, wtsum, cpu0, cpu
  logical :: overlap, use_weights

  namelist /nml/ temperature, mu, nmax, nstep, dr_max

  write(output_unit,'(a)') 'mc_muca'
  write(output_unit,'(a)') 'Grand canonical / multicanonical simulation of system of atoms'
  call introduction

  call random_init ( .false., .true. ) ! Initialize random number generator

  ! Set sensible default run parameters for testing
  temperature = 0.95
  mu          = -3.14
  nmax        = 285
  nstep       = 10000000
  dr_max      = 0.15

  ! Read run parameters from namelist
  ! Comment out, or replace, this section if you don't like namelists
  read ( unit=input_unit, nml=nml, iostat=ioerr )
  if ( ioerr /= 0 ) then
     write ( unit=error_unit, fmt='(a,i10)') 'Error reading namelist nml from standard input', ioerr
     if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
     if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
     stop 'Error in muca'
  end if

  beta = 1.0 / temperature

  ! Write out run parameters
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Number of steps        ',  nstep
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Maximum number of atoms',  nmax
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Temperature T',            temperature
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Maximum displacement',     dr_max
  write ( unit=output_unit, fmt='(a,t31,f10.4)' ) 'Chemical potential mu',    mu

  ! Allocate arrays
  allocate ( r(3,nmax) )
  allocate ( h(0:nmax), phi(0:nmax) )
  allocate ( mu_old(nmax), wt_old(nmax), mu_new(nmax), wt_new(nmax) )

  ! Read mu table and statistical weights, if the old HDF5 file exists
  inquire ( file=hdf5_old, exist=use_weights )

  if ( use_weights ) then

     write ( unit=output_unit, fmt='(a)' ) 'Multicanonical simulation'

     ! Open existing HDF5 file, for reading only
     write ( unit=output_unit, fmt='(a,t31,a)' ) 'Opening HDF5 file', hdf5_old
     call open_file ( hdf5_old, readonly=.true. )

     ! Read nmax from attributes and confirm that it matches, for simplicity
     call read_attribute ( 'Nmax', n )
     if ( n /= nmax ) stop 'We insist on a matching value of Nmax in this file'

     ! Read mu table and statistical weights from datasets
     write ( unit=output_unit, fmt='(a)' ) 'Reading multicanonical mu and weights '
     call read_dataset ( 'mu', mu_old )
     call read_dataset ( 'wt', wt_old )

     ! Close HDF5 file
     call close_file

  else

     ! Set mu table for standard grand canonical simulation
     ! These values will persist through successive simulations until over-written 
     write ( unit=output_unit, fmt='(a)' ) 'Grand canonical simulation'
     mu_old = mu
     wt_old = 0
  end if

  ! phi(N) = F(N) where F(N) = F(N-1) + mu(N)
  ! Uses cumulative sum function defined below
  phi(0:nmax)  = [ 0.0, cumsum(mu_old) ]

  ! Read initial configuration
  r = 0. ! Fill with zeros
  n = n_config ( config_old )                   ! Get n (can be zero)
  call read_config ( config_old, box, r(:,:n) ) ! Get box and r
  write ( unit=output_unit, fmt='(a,t31,i10)'   ) 'Initial n',   n
  write ( unit=output_unit, fmt='(a,t31,3f10.4)') 'Initial box', box

  ! Check for overlaps; calculate total potential energy & virial
  call pot_all ( box, r(:,:n), u, w, overlap )
  if ( overlap ) stop 'Overlaps detected in initial configuration'

  ! Initialize histograms and counters
  call cpu_time ( cpu0 )
  h = 0.
  r_accepted = 0
  n_accepted = 0

  ! Column headings
  write ( unit=output_unit, fmt='(4a10)' ) 'Step', 'CPU', 'Move %', 'N move %' 

  do step = 1, nstep ! Begin loop over steps

     call r_move ( beta, dr_max, box, r(:,:n), u, w, r_accepted )

     call n_move ( beta, phi, box, n, r, u, w, n_accepted )

     ! Accumulate histogram
     if ( n <= nmax ) h(n) = h(n) + 1.

     ! Write out move acceptance rates
     if ( trigger(step) ) then
        r_percent = 100*real(r_accepted)/real(step)
        n_percent = 100*real(n_accepted)/real(step)
        call cpu_time ( cpu )
        write ( unit=output_unit, fmt='(i10,3f10.2)' ) step, cpu-cpu0, r_percent, n_percent
     end if

  end do ! End loop over steps

  call cpu_time ( cpu )

  ! Write final configuration
  call write_config ( config_new, box, r(:,:n) )

  ! Calculate new estimates of mu and statistical weights
  do n = 1, nmax
     mu_new(n) = 0
     wt_new(n) = 0
     if ( h(n-1)*h(n) > tiny(1.0) ) then
        mu_new(n) = -temperature*log(h(n)/h(n-1)) + phi(n) - phi(n-1)
        wt_new(n) = h(n-1)*h(n)/(h(n-1)+h(n))
     end if
  end do

  ! Mix new and old mu and statistical weights
  do n = 1, nmax
     wtsum = wt_old(n) + wt_new(n)
     if ( wtsum > tiny(1.0) ) then
        mu_old(n) = (wt_old(n)*mu_old(n) + wt_new(n)*mu_new(n)) / wtsum
        wt_old(n) = wtsum
     end if
  end do
  mu_new = mu_old
  wt_new = wt_old

  ! Create new HDF5 file, write simulation parameters as attributes
  call open_file ( hdf5_new )
  call write_attribute ( 'Title', 'Monte Carlo MUCA' )
  call write_attribute ( 'Nmax', nmax )
  call write_attribute ( 'nstep', nstep )
  call write_attribute ( 'L', box )
  call write_attribute ( 'V', product(box) )
  call write_attribute ( 'T', temperature )
  call write_attribute ( 'mu', mu )

  ! Write histogram, mu table, and statistical weights as datasets
  call write_dataset   ( 'h', h )
  call write_dataset   ( 'mu', mu_new )
  call write_dataset   ( 'wt', wt_new )

  ! Close HDF5 file
  call close_file

  deallocate ( r )
  deallocate ( h, phi )
  deallocate ( mu_old, wt_old, mu_new, wt_new )
  write ( unit=output_unit, fmt='(a,f10.2)') 'Program ends, CPU = ', cpu-cpu0

contains

  function cumsum ( x ) result ( c ) ! Cumulative sum function for supplied array
    real, dimension(:),      intent(in) :: x
    real, dimension(size(x))            :: c

    integer :: j

    if ( size(x) == 0 ) return

    c(1) = x(1)
    do j = 2, size(x)
       c(j) = c(j-1) + x(j)
    end do

  end function cumsum

end program mc_muca

