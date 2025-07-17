! hdf5_module.f90
! Routines for writing simulation results to HDF5 file
! Compile with -fdefault-real-8 or similar, for 8-byte default real variables
module hdf5_module

  ! DISCLAIMER
  ! (c) 2022 M P Allen, exclusively for the CCP5 Summer School, for educational purposes only.
  ! The author makes no warranties about, and disclaims liability for all uses of, this software.
  ! The author does not recommend use of this software for any purpose.

  ! High-level routines to write out simulation results to an HDF5 file for further analysis.

  ! We make a number of assumptions and restrictions, so as to keep the interface simple.
  ! We only use a flat HDF5 file: no groups.
  ! Only one file is in use at any one time: it is opened, written to or read from, and closed.
  ! This allows the file_id to be private to this module.
  ! The file contains attributes and datasets.
  ! The attributes are attached to the root of the file, and contain simulation parameters.
  ! The datasets typically contain step-by-step simulation values.
  ! In both cases we only consider a limited subset of possible ranks, dimensions, and data kinds.
  ! All dimensions etc are determined by the variables passed to the routines.
  ! Be aware that HDF5 is much more flexible than this!

  ! We do not test the error flag returned by each HDF5 routine, to keep things simple.
  ! Obviously, in a practical application, better practice is to test these error flags.

  use HDF5
  implicit none
  private
  public :: open_file, close_file, read_dataset, read_attribute, write_dataset, write_attribute

  ! File identifier and file_open flag. These are saved as module variables.
  integer(HID_T) :: file_id
  logical        :: file_open = .false.

  ! Generic procedures allowing a simple user interface

  interface write_dataset
     module procedure write_dataset_integer_1
     module procedure write_dataset_integer_2
     module procedure write_dataset_real_1
     module procedure write_dataset_real_2
     module procedure write_dataset_real_3
  end interface write_dataset

  interface read_dataset
     module procedure read_dataset_integer_1
     module procedure read_dataset_real_1
  end interface read_dataset

  interface write_attribute
     module procedure write_attribute_integer_0
     module procedure write_attribute_integer_1
     module procedure write_attribute_real_0
     module procedure write_attribute_real_1
     module procedure write_attribute_string
  end interface write_attribute

  interface read_attribute
     module procedure read_attribute_integer_0
     module procedure read_attribute_integer_1
     module procedure read_attribute_real_0
     module procedure read_attribute_real_1
     module procedure read_attribute_string
  end interface read_attribute

contains

  subroutine open_file ( filename, readonly )
    implicit none
    character(len=*),  intent(in) :: filename ! Supplied file name
    logical, optional, intent(in) :: readonly ! Optionally indicates read only access

    integer        :: error     ! HDF5 error flag
    integer(HID_T) :: f_id      ! File id returned by HDF5 routine
    logical        :: read_only ! Local copy of readonly argument

    if ( file_open ) stop 'Attempt to open 2 HDF5 files at once'
    
    call h5open_f ( error ) ! Open Fortran interface

    read_only = .false.
    if ( present(readonly) ) read_only = readonly

    if ( read_only ) then
       call h5fopen_f ( filename, H5F_ACC_RDONLY_F, f_id, error ) ! Open existing file, define file id
    else
       call h5fcreate_f ( filename, H5F_ACC_TRUNC_F, f_id, error ) ! Create new file, define file id
    end if
    file_id   = f_id
    file_open = .true.

  end subroutine open_file

  subroutine close_file
    implicit none

    integer :: error ! HDF5 error flag

    if ( .not. file_open ) stop 'Attempt to close an HDF5 file that is not open'
    
    call h5fclose_f ( file_id, error ) ! Close file
    call h5close_f  ( error ) ! Close Fortran interface
    file_open = .false.
    
  end subroutine close_file

  subroutine read_dataset_integer_1 ( name, dataset )
    implicit none
    character(len=*),               intent(in)  :: name    ! Dataset name
    integer,          dimension(:), intent(out) :: dataset ! Dataset values

    integer(HID_T)                 :: dataset_id   ! Dataset identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Dataset dimensions
    integer                        :: rank         ! Dataset rank
    integer                        :: error        ! Error flag

    dims = shape ( dataset )
    rank = 1
    call h5dopen_f  ( file_id, name, dataset_id, error )
    call h5dread_f  ( dataset_id, H5T_NATIVE_INTEGER, dataset, dims, error )
    call h5dclose_f ( dataset_id, error )

  end subroutine read_dataset_integer_1

  subroutine write_dataset_integer_1 ( name, dataset )
    implicit none
    character(len=*),               intent(in) :: name    ! Dataset name
    integer,          dimension(:), intent(in) :: dataset ! Dataset values

    integer(HID_T)                 :: dataset_id   ! Dataset identifier
    integer(HID_T)                 :: dataspace_id ! Dataspace identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Dataset dimensions
    integer                        :: rank         ! Dataset rank
    integer                        :: error        ! Error flag

    dims = shape ( dataset )
    rank = 1
    call h5screate_simple_f ( rank, dims, dataspace_id, error )
    call h5dcreate_f ( file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error )
    call h5dwrite_f  ( dataset_id, H5T_NATIVE_INTEGER, dataset, dims, error )
    call h5dclose_f  ( dataset_id, error )
    call h5sclose_f  ( dataspace_id, error )

  end subroutine write_dataset_integer_1

  subroutine write_dataset_integer_2 ( name, dataset )
    implicit none
    character(len=*),                 intent(in) :: name    ! Dataset name
    integer,          dimension(:,:), intent(in) :: dataset ! Dataset values

    integer(HID_T)                 :: dataset_id   ! Dataset identifier
    integer(HID_T)                 :: dataspace_id ! Dataspace identifier
    integer(HSIZE_T), dimension(2) :: dims         ! Dataset dimensions
    integer                        :: rank         ! Dataset rank
    integer                        :: error        ! Error flag

    dims = shape ( dataset )
    rank = 2
    call h5screate_simple_f ( rank, dims, dataspace_id, error )
    call h5dcreate_f ( file_id, name, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error )
    call h5dwrite_f  ( dataset_id, H5T_NATIVE_INTEGER, dataset, dims, error )
    call h5dclose_f  ( dataset_id, error )
    call h5sclose_f  ( dataspace_id, error )

  end subroutine write_dataset_integer_2

  subroutine read_dataset_real_1 ( name, dataset )
    implicit none
    character(len=*),               intent(in)  :: name    ! Dataset name
    real,             dimension(:), intent(out) :: dataset ! Dataset values

    integer(HID_T)                 :: dataset_id   ! Dataset identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Dataset dimensions
    integer                        :: rank         ! Dataset rank
    integer                        :: error        ! Error flag

    ! NB we assume the compiler option -fdefault-real-8 or similar
    ! so that 8-byte default real variables map onto H5T_NATIVE_DOUBLE
    
    dims = shape ( dataset )
    rank = 1
    call h5dopen_f  ( file_id, name, dataset_id, error )
    call h5dread_f  ( dataset_id, H5T_NATIVE_DOUBLE, dataset, dims, error )
    call h5dclose_f ( dataset_id, error )

  end subroutine read_dataset_real_1

  subroutine write_dataset_real_1 ( name, dataset )
    implicit none
    character(len=*),               intent(in) :: name    ! Dataset name
    real,             dimension(:), intent(in) :: dataset ! Dataset values

    integer(HID_T)                 :: dataset_id   ! Dataset identifier
    integer(HID_T)                 :: dataspace_id ! Dataspace identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Dataset dimensions
    integer                        :: rank         ! Dataset rank
    integer                        :: error        ! Error flag

    ! NB we assume the compiler option -fdefault-real-8 or similar
    ! so that 8-byte default real variables map onto H5T_NATIVE_DOUBLE
    
    dims = shape ( dataset )
    rank = 1
    call h5screate_simple_f ( rank, dims, dataspace_id, error )
    call h5dcreate_f ( file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error )
    call h5dwrite_f  ( dataset_id, H5T_NATIVE_DOUBLE, dataset, dims, error )
    call h5dclose_f  ( dataset_id, error )
    call h5sclose_f  ( dataspace_id, error )

  end subroutine write_dataset_real_1

  subroutine write_dataset_real_2 ( name, dataset )
    implicit none
    character(len=*),                 intent(in) :: name    ! Dataset name
    real,             dimension(:,:), intent(in) :: dataset ! Dataset values

    integer(HID_T)                 :: dataset_id   ! Dataset identifier
    integer(HID_T)                 :: dataspace_id ! Dataspace identifier
    integer(HSIZE_T), dimension(2) :: dims         ! Dataset dimensions
    integer                        :: rank         ! Dataset rank
    integer                        :: error        ! Error flag

    ! NB we assume the compiler option -fdefault-real-8 or similar
    ! so that 8-byte default real variables map onto H5T_NATIVE_DOUBLE

    dims = shape ( dataset )
    rank = 2
    call h5screate_simple_f ( rank, dims, dataspace_id, error )
    call h5dcreate_f ( file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error )
    call h5dwrite_f  ( dataset_id, H5T_NATIVE_DOUBLE, dataset, dims, error )
    call h5dclose_f  ( dataset_id, error )
    call h5sclose_f  ( dataspace_id, error )

  end subroutine write_dataset_real_2

  subroutine write_dataset_real_3 ( name, dataset )
    implicit none
    character(len=*),                   intent(in) :: name    ! Dataset name
    real,             dimension(:,:,:), intent(in) :: dataset ! Dataset values

    integer(HID_T)                 :: dataset_id   ! Dataset identifier
    integer(HID_T)                 :: dataspace_id ! Dataspace identifier
    integer(HSIZE_T), dimension(3) :: dims         ! Dataset dimensions
    integer                        :: rank         ! Dataset rank
    integer                        :: error        ! Error flag

    ! NB we assume the compiler option -fdefault-real-8 or similar
    ! so that 8-byte default real variables map onto H5T_NATIVE_DOUBLE

    dims = shape ( dataset )
    rank = 3
    call h5screate_simple_f ( rank, dims, dataspace_id, error )
    call h5dcreate_f ( file_id, name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error )
    call h5dwrite_f  ( dataset_id, H5T_NATIVE_DOUBLE, dataset, dims, error )
    call h5dclose_f  ( dataset_id, error )
    call h5sclose_f  ( dataspace_id, error )

  end subroutine write_dataset_real_3

  subroutine read_attribute_integer_0 ( name, attribute )
    implicit none
    character(len=*), intent(in)  :: name      ! Attribute name
    integer,          intent(out) :: attribute ! Attribute value

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: error        ! Error flag

    call h5aopen_f  ( file_id, name, attribute_id, error )
    call h5aread_f  ( attribute_id, H5T_NATIVE_INTEGER, attribute, dims, error )
    call h5aclose_f ( attribute_id, error )

  end subroutine read_attribute_integer_0

  subroutine write_attribute_integer_0 ( name, attribute )
    implicit none
    character(len=*), intent(in) :: name      ! Attribute name
    integer,          intent(in) :: attribute ! Attribute value

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HID_T)                 :: attrspace_id ! Attribute-space identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: error        ! Error flag

    dims = [ 0 ]
    call h5screate_f ( H5S_SCALAR_F, attrspace_id, error )
    call h5acreate_f ( file_id, name, H5T_NATIVE_INTEGER, attrspace_id, attribute_id, error )
    call h5awrite_f  ( attribute_id, H5T_NATIVE_INTEGER, attribute, dims, error )
    call h5aclose_f  ( attribute_id, error )
    call h5sclose_f  ( attrspace_id, error )

  end subroutine write_attribute_integer_0

  subroutine read_attribute_integer_1 ( name, attribute )
    implicit none
    character(len=*),               intent(in)  :: name      ! Attribute name
    integer,          dimension(:), intent(out) :: attribute ! Attribute values

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: rank         ! Attribute rank
    integer                        :: error        ! Error flag

    rank = 1
    dims = shape ( attribute )
    call h5aopen_f  ( file_id, name, attribute_id, error )
    call h5aread_f  ( attribute_id, H5T_NATIVE_INTEGER, attribute, dims, error )
    call h5aclose_f ( attribute_id, error )

  end subroutine read_attribute_integer_1

  subroutine write_attribute_integer_1 ( name, attribute )
    implicit none
    character(len=*),               intent(in) :: name      ! Attribute name
    integer,          dimension(:), intent(in) :: attribute ! Attribute values

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HID_T)                 :: attrspace_id ! Attribute-space identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: rank         ! Attribute rank
    integer                        :: error        ! Error flag

    rank = 1
    dims = shape ( attribute )
    call h5screate_simple_f ( rank, dims, attrspace_id, error )
    call h5acreate_f ( file_id, name, H5T_NATIVE_INTEGER, attrspace_id, attribute_id, error )
    call h5awrite_f  ( attribute_id, H5T_NATIVE_INTEGER, attribute, dims, error )
    call h5aclose_f  ( attribute_id, error )
    call h5sclose_f  ( attrspace_id, error )

  end subroutine write_attribute_integer_1

  subroutine read_attribute_real_0 ( name, attribute )
    implicit none
    character(len=*), intent(in)  :: name      ! Attribute name
    real,             intent(out) :: attribute ! Attribute value

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: error        ! Error flag

    ! NB we assume the compiler option -fdefault-real-8 or similar
    ! so that 8-byte default real variables map onto H5T_NATIVE_DOUBLE

    call h5aopen_f  ( file_id, name, attribute_id, error )
    call h5aread_f  ( attribute_id, H5T_NATIVE_DOUBLE, attribute, dims, error )
    call h5aclose_f ( attribute_id, error )

  end subroutine read_attribute_real_0

  subroutine write_attribute_real_0 ( name, attribute )
    implicit none
    character(len=*), intent(in) :: name      ! Attribute name
    real,             intent(in) :: attribute ! Attribute value

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HID_T)                 :: attrspace_id ! Attribute-space identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: error        ! Error flag

    ! NB we assume the compiler option -fdefault-real-8 or similar
    ! so that 8-byte default real variables map onto H5T_NATIVE_DOUBLE

    dims = [ 0 ]
    call h5screate_f ( H5S_SCALAR_F, attrspace_id, error )
    call h5acreate_f ( file_id, name, H5T_NATIVE_DOUBLE, attrspace_id, attribute_id, error )
    call h5awrite_f  ( attribute_id, H5T_NATIVE_DOUBLE, attribute, dims, error )
    call h5aclose_f  ( attribute_id, error )
    call h5sclose_f  ( attrspace_id, error )

  end subroutine write_attribute_real_0

  subroutine read_attribute_real_1 ( name, attribute )
    implicit none
    character(len=*),               intent(in)  :: name      ! Attribute name
    real,             dimension(:), intent(out) :: attribute ! Attribute values

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: rank         ! Attribute rank
    integer                        :: error        ! Error flag

    ! NB we assume the compiler option -fdefault-real-8 or similar
    ! so that 8-byte default real variables map onto H5T_NATIVE_DOUBLE

    rank = 1
    dims = shape ( attribute )
    call h5aopen_f  ( file_id, name, attribute_id, error )
    call h5aread_f  ( attribute_id, H5T_NATIVE_DOUBLE, attribute, dims, error )
    call h5aclose_f ( attribute_id, error )

  end subroutine read_attribute_real_1

  subroutine write_attribute_real_1 ( name, attribute )
    implicit none
    character(len=*),               intent(in) :: name      ! Attribute name
    real,             dimension(:), intent(in) :: attribute ! Attribute values

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HID_T)                 :: attrspace_id ! Attribute-space identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: rank         ! Attribute rank
    integer                        :: error        ! Error flag

    ! NB we assume the compiler option -fdefault-real-8 or similar
    ! so that 8-byte default real variables map onto H5T_NATIVE_DOUBLE

    rank = 1
    dims = shape ( attribute )
    call h5screate_simple_f ( rank, dims, attrspace_id, error )
    call h5acreate_f ( file_id, name, H5T_NATIVE_DOUBLE, attrspace_id, attribute_id, error )
    call h5awrite_f  ( attribute_id, H5T_NATIVE_DOUBLE, attribute, dims, error )
    call h5aclose_f  ( attribute_id, error )
    call h5sclose_f  ( attrspace_id, error )

  end subroutine write_attribute_real_1

  subroutine read_attribute_string ( name, attribute )
    implicit none
    character(len=*), intent(in)  :: name      ! Attribute name
    character(len=*), intent(out) :: attribute ! Attribute value

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HID_T)                 :: chartype_id  ! Character type identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: error        ! Error flag

    dims = [ len(attribute, kind=HID_T) ]
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, chartype_id, error )
    call h5tset_size_f ( chartype_id, dims(1), error )
    call h5aopen_f ( file_id, name, attribute_id, error )
    call h5aread_f  ( attribute_id, chartype_id, attribute, dims, error )
    call h5aclose_f  ( attribute_id, error )

  end subroutine read_attribute_string

  subroutine write_attribute_string ( name, attribute )
    implicit none
    character(len=*), intent(in) :: name      ! Attribute name
    character(len=*), intent(in) :: attribute ! Attribute value

    integer(HID_T)                 :: attribute_id ! Attribute identifier
    integer(HID_T)                 :: attrspace_id ! Attribute-space identifier
    integer(HID_T)                 :: chartype_id  ! Character type identifier
    integer(HSIZE_T), dimension(1) :: dims         ! Attribute dimensions
    integer                        :: error        ! Error flag

    dims = [ len(attribute, kind=HID_T) ]
    call h5tcopy_f ( H5T_NATIVE_CHARACTER, chartype_id, error )
    call h5tset_size_f ( chartype_id, dims(1), error )
    call h5screate_f ( H5S_SCALAR_F, attrspace_id, error )
    call h5acreate_f ( file_id, name, chartype_id, attrspace_id, attribute_id, error )
    call h5awrite_f  ( attribute_id, chartype_id, attribute, dims, error )
    call h5aclose_f  ( attribute_id, error )
    call h5sclose_f  ( attrspace_id, error )

  end subroutine write_attribute_string

end module hdf5_module
