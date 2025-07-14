! config_io_module.f90
! Routines for atomic configuration input and output
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module config_io_module

  ! DISCLAIMER
  ! (c) 2022, 2024 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  use, intrinsic :: iso_fortran_env, only : error_unit, iostat_end, iostat_eor

  implicit none
  private

  ! Public routines
  public :: n_config, read_config, write_config

contains

  function n_config ( file ) result ( n )
    implicit none
    character(len=*), intent(in) :: file ! Supplied filename
    integer                      :: n    ! Number of atoms

    ! Reads number of atoms from file, which must exist
    ! This allows us to allocate the appropriate arrays in the main program
    ! For us this is more convenient than allocating them in this module

    integer :: unit, ioerr

    ! Open given file, will terminate on any errors

    open ( newunit=unit, file=file, status='old', action='read', iostat=ioerr )
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error opening ', file, ioerr
       stop 'Error in n_config'
    end if

    ! Read number of atoms from first record

    read ( unit=unit, fmt=*, iostat=ioerr ) n
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading n from ', file, ioerr
       if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
       if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
       stop 'Error in n_config'
    end if

    ! n must be positive or zero
    if ( n < 0 ) stop 'n error in n_config'

    close ( unit=unit )

  end function n_config

  subroutine read_config ( file, box, r, p )
    implicit none
    character(len=*),               intent(in)  :: file ! Supplied filename
    real, dimension(3),             intent(out) :: box  ! Simulation box lengths
    real, dimension(:,:),           intent(out) :: r    ! Atomic positions (3,n)
    real, dimension(:,:), optional, intent(out) :: p    ! Atomic momenta (3,n)

    ! Reads atomic configuration from file, which must exist
    ! It is possible for n to be zero, in which case no coordinates will be read

    integer :: unit, ioerr, i, n

    ! Open given file, will terminate on any errors

    open ( newunit=unit, file=file, status='old', action='read', iostat=ioerr )
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error opening ', file, ioerr
       stop 'Error in read_config'
    end if

    ! Read number of atoms from first record

    read ( unit=unit, fmt=*, iostat=ioerr ) n
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading n from ', file, ioerr
       if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
       if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
       stop 'Error in read_config'
    end if

    if ( n < 0 ) stop 'n error in read_config'
    if ( any ( shape(r) /= [3,n] ) ) stop 'r dimension error in read_config'
    if ( present(p) ) then
       if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in write_config'
    end if

    ! Read box lengths from second record

    read ( unit=unit, fmt=*, iostat=ioerr ) box(:)
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading box(:) from ', file, ioerr
       if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
       if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
       stop 'Error in read_config'
    end if
    
    ! If n==0 there are no coordinates to read
    if ( n == 0 ) then
       close ( unit = unit )
       return
    end if

    ! Expected format for the rest of the file is one record per atom
    ! containing either r(:,i), p(:,i) for MD, or just r(:,i) for MC

    if ( present(p) ) then

       do i = 1, n
          read ( unit=unit, fmt=*, iostat=ioerr ) r(:,i), p(:,i)
          if ( ioerr /= 0 ) then
             write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading r, p from ', file, ioerr
             if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
             if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
             stop 'Error in read_config'
          end if
       end do

    else

       do i = 1, n
          read ( unit=unit, fmt=*, iostat=ioerr ) r(:,i)
          if ( ioerr /= 0 ) then
             write ( unit=error_unit, fmt='(a,a,i10)') 'Error reading r from ', file, ioerr
             if ( ioerr == iostat_eor ) write ( unit=error_unit, fmt='(a)') 'End of record'
             if ( ioerr == iostat_end ) write ( unit=error_unit, fmt='(a)') 'End of file'
             stop 'Error in read_config'
          end if
       end do

    end if

    close ( unit=unit )

  end subroutine read_config

  subroutine write_config ( file, box, r, p )
    implicit none
    character(len=*),               intent(in) :: file ! Supplied filename
    real, dimension(3),             intent(in) :: box  ! Simulation box lengths
    real, dimension(:,:),           intent(in) :: r    ! Atomic positions (3,n)
    real, dimension(:,:), optional, intent(in) :: p    ! Atomic momenta (3,n)

    ! Writes out atomic configuration to file, which will be overwritten if it exists
    
    integer :: unit, ioerr, i, n

    if ( size(r,dim=1) /= 3 ) stop 'r dimension error in write_config'
    n =  size(r,dim=2)
    if ( present(p) ) then
       if ( any ( shape(p) /= [3,n] ) ) stop 'p dimension error in write_config'
    end if

    ! Open given file, replacing it if it already exists, will terminate on any errors
    open ( newunit=unit, file=file, status='replace', iostat=ioerr )
    if ( ioerr /= 0 ) then
       write ( unit=error_unit, fmt='(a,a,i10)') 'Error opening ', file, ioerr
       stop 'Error in write_config'
    end if

    ! Write number of atoms to first record, box length to second record
    write ( unit=unit, fmt='(i15)'    ) n
    write ( unit=unit, fmt='(3f15.10)') box(:)

    if ( present ( p ) ) then
       do i = 1, n
          write ( unit=unit, fmt='(*(f15.10))' ) r(:,i), p(:,i)
       end do
    else
       do i = 1, n
          write ( unit=unit, fmt='(*(f15.10))' ) r(:,i)
       end do
    end if

    close ( unit=unit )

  end subroutine write_config

end module config_io_module
