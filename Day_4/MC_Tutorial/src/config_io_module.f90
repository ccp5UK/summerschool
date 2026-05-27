! config_io_module.f90
! Routines for atomic MC configuration input and output
module config_io_module

  ! DISCLAIMER
  ! (c) 2022, 2026 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, iostat_end, iostat_eor, dp => real64

  implicit none

  ! Public routines
  private
  public :: read_config, write_config

  ! Generic interface for MC or MD functions
  interface read_config
     module procedure read_config_mc ! pbc box and positions
     module procedure read_config_md ! positions and velocities
  end interface read_config
  interface write_config
     module procedure write_config_mc ! pbc box and positions
     module procedure write_config_md ! positions and velocities
  end interface write_config
contains

  subroutine read_config_mc ( file, box, r )
    implicit none
    character(*),          intent(in)  :: file   ! Supplied filename
    real(dp),              intent(out) :: box(3) ! Simulation box lengths
    real(dp), allocatable, intent(out) :: r(:,:) ! Atomic positions (3,n)

    ! Reads atomic MC configuration (box dimensions and coordinates) from file, which must exist
    ! Allocates atomic position array, r(3,n), where n is the number of atoms
    ! It is possible for n to be zero, in which case no coordinates will be read

    ! The file must be in a very specific version of extended xyz (extxyz) format
    ! See https://github.com/libAtoms/extxyz
    ! We are not intending to parse a general extxyz file here

    integer            :: unit, ioerr, i, j, n
    real(dp)           :: cell(3,3)
    character(len=200) :: line
    character(len=1)   :: s

    ! Open given file, will terminate on any errors

    write ( unit=output_unit, fmt='(a,t31,a)') 'Reading MC config from file', file
    open ( newunit=unit, file=file, status='old', action='read', iostat=ioerr )
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error opening ', file, ioerr
       stop 'Error in read_config_mc'
    end if

    ! Read number of atoms from first record

    read ( unit=unit, fmt=*, iostat=ioerr ) n
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading n from ', file, ioerr
       if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
       if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
       stop 'Error in read_config_mc'
    end if

    ! Allocate atomic position array

    if ( n < 0 ) stop 'n error in read_config_mc'
    allocate ( r(3,n) )

    ! Read box lengths from second record
    ! Extended xyz format has a 9-element (3x3) Lattice representing the general orthorhombic cell
    ! but we expect this to be diagonal, for a cuboidal cell with 3 box lengths
    ! We delimit the Lattice elements with [...] which is one of the allowed extxyz options

    read ( unit=unit, fmt='(a)', iostat=ioerr ) line
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading line from ', file, ioerr
       if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
       if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
       stop 'Error in read_config_mc'
    end if
    if ( index(line,'Lattice=[') == 0 ) then
       write ( unit=error_unit, fmt='(a)') 'Expected Lattice=[...] for box dimensions'
       stop 'Error in read_config_mc'
    end if
    i = index(line,'[')
    j = index(line,']')
    if ( i==0 .or. j==0 ) then
       write ( unit=error_unit, fmt='(a)') 'Expected [...] for box dimensions'
       stop 'Error in read_config_mc'
    end if
    read( unit=line(i+1:j-1), fmt=*, iostat=ioerr) cell
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Expected 9 reals in ', line(i+1:j-1), ioerr
       if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
       if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
       stop 'Error in read_config_mc'
    end if
    do i = 1, 3
       do j = 1, 3
          if ( i == j ) then
             box(i) = cell(i,i)
          else
             if ( abs(cell(i,j))>tiny(1.0_dp) )  then
                write ( unit=error_unit, fmt='(a,e15.6)') 'Expected cell to be diagonal', cell(i,j)
                stop 'Error in read_config_mc'
             end if
          end if
       end do
    end do

    ! Expected format for the rest of the file is one record per atom containing s(i) r(:,i) for MC

    do i = 1, n
       read ( unit=unit, fmt=*, iostat=ioerr ) s, r(:,i)
       if ( ioerr /= 0 ) then
          write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading r from ', file, ioerr
          if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
          if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
          stop 'Error in read_config_mc'
       end if
    end do

    close ( unit=unit )

  end subroutine read_config_mc

  subroutine write_config_mc ( file, box, r )
    implicit none
    character(*), intent(in) :: file   ! Supplied filename
    real(dp),     intent(in) :: box(3) ! Simulation box lengths
    real(dp),     intent(in) :: r(:,:) ! Atomic positions (3,n)

    ! Writes out atomic configuration to file, which will be overwritten if it exists
    ! The file conforms with extended xyz format, see https://github.com/libAtoms/extxyz

    integer  :: unit, ioerr, i, n
    real(dp) :: cell(3,3)

    character(len=*), parameter :: head = 'Lattice=['
    character(len=*), parameter :: tail = '] Properties=species:S:1:pos:R:3'
    character(len=*), parameter :: form = '(a,f15.10,3f4.1,f15.10,3f4.1,f15.10,a)'

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in write_config_mc'
    n =  size(r,dim=2)

    ! Open given file, replacing it if it already exists, will terminate on any errors

    write ( unit=output_unit, fmt='(a,t31,a)') 'Writing MC config to file',file
    open ( newunit=unit, file=file, status='replace', iostat=ioerr )
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error opening ', file, ioerr
       stop 'Error in write_config_mc'
    end if

    ! Write number of atoms to first record

    write ( unit=unit, fmt='(i15)' ) n

    ! Write box length to second record (diagonal 3x3 cell)

    cell = 0.0_dp
    do i = 1, 3
       cell(i,i) = box(i)
    end do
    write ( unit=unit, fmt=form ) head, cell, tail

    ! One record per atom containing s(i) r(:,i) for MC

    do i = 1, n
       write ( unit=unit, fmt='(a1,3f15.10)' ) 'X', r(:,i)
    end do

    close ( unit=unit )

  end subroutine write_config_mc

  subroutine read_config_md ( file, r, p )
    implicit none
    character(*),          intent(in)  :: file   ! Supplied filename
    real(dp), allocatable, intent(out) :: r(:,:) ! Atomic positions (3,n)
    real(dp), allocatable, intent(out) :: p(:,:) ! Atomic momenta (3,n)

    ! Reads atomic MD configuration (coordinates and momenta) from file, which must exist
    ! Allocates atomic position and momentum arrays, r(3,n), p(3,n), where n is the number of atoms
    ! It is possible for n to be zero, in which case no coordinates will be read
    ! There is no periodic box

    ! The file must be in a very specific version of extended xyz (extxyz) format
    ! See https://github.com/libAtoms/extxyz
    ! We are not intending to parse a general extxyz file here

    integer            :: unit, ioerr, i, n
    character(len=200) :: line
    character(len=1)   :: s

    ! Open given file, will terminate on any errors

    write ( unit=output_unit, fmt='(a,t31,a)') 'Reading MD config from file', file
    open ( newunit=unit, file=file, status='old', action='read', iostat=ioerr )
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error opening ', file, ioerr
       stop 'Error in read_config_md'
    end if

    ! Read number of atoms from first record

    read ( unit=unit, fmt=*, iostat=ioerr ) n
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading n from ', file, ioerr
       if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
       if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
       stop 'Error in read_config_md'
    end if

    ! Allocate atomic position and momentum arrays

    if ( n < 0 ) stop 'n error in read_config_md'
    allocate ( r(3,n), p(3,n) )

    ! Second record contains no information for us, but check Lattice/cell/box is not present

    read ( unit=unit, fmt='(a)', iostat=ioerr ) line
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading line from ', file, ioerr
       if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
       if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
       stop 'Error in read_config_md'
    end if
    if ( index(line,'Lattice=[') /= 0 ) then
       write ( unit=error_unit, fmt='(a)') 'Did NOT expect Lattice=[...] for box dimensions'
       stop 'Error in read_config_md'
    end if

    ! Expected format for the rest of the file is one record per atom containing s(i) r(:,i) p(:,i) for MD

    do i = 1, n
       read ( unit=unit, fmt=*, iostat=ioerr ) s, r(:,i), p(:,i)
       if ( ioerr /= 0 ) then
          write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading r, p from ', file, ioerr
          if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
          if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
          stop 'Error in read_config_md'
       end if
    end do

    close ( unit=unit )

  end subroutine read_config_md

  subroutine write_config_md ( file, r, p )
    implicit none
    character(*), intent(in) :: file   ! Supplied filename
    real(dp),     intent(in) :: r(:,:) ! Atomic positions (3,n)
    real(dp),     intent(in) :: p(:,:) ! Atomic momenta (3,n)

    ! Writes out atomic configuration to file, which will be overwritten if it exists
    ! The file conforms with extended xyz format, see https://github.com/libAtoms/extxyz
    ! There is no periodic box

    integer :: unit, ioerr, i, n

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in write_config_md'
    n =  size(r,dim=2)
    if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in write_config_md'

    ! Open given file, replacing it if it already exists, will terminate on any errors

    write ( unit=output_unit, fmt='(a,t31,a)') 'Writing MD config to file',file
    open ( newunit=unit, file=file, status='replace', iostat=ioerr )
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error opening ', file, ioerr
       stop 'Error in write_config_md'
    end if

    ! Write number of atoms to first record

    write ( unit=unit, fmt='(i15)' ) n

    ! Write standard properties string to second record

    write ( unit=unit, fmt='(a)' ) 'Properties=species:S:1:pos:R:3:velo:R:3'

    ! One record per atom containing s(i) r(:,i) p(:,i) for MD

    do i = 1, n
       write ( unit=unit, fmt='(a1,6f15.10)' ) 'X', r(:,i), p(:,i)
    end do

    close ( unit=unit )

  end subroutine write_config_md

end module config_io_module
