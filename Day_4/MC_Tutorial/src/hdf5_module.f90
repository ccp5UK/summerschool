! hdf5_module_2.f90
! Routines for writing simulation results to HDF5 file

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
  ! In both cases we only consider a limited subset of possible  data kinds.
  ! We only use character strings in the attributes, not the datasets,
  ! and they are scalars (of arbitrary length), not arrays.
  ! All ranks and dimensions etc are determined by the variables passed to the routines.
  ! Hence, on reading data, the dataset names, ranks and dimensions must be known by the calling program.
  ! Be aware that HDF5 is much more flexible than this!

  ! We do not test the error flag returned by each HDF5 routine, to keep things simple.
  ! Obviously, in a practical application, better practice is to test these error flags.

  ! Working precision for reals, dp, is real64 throughout
  use, intrinsic :: iso_fortran_env, only : output_unit, dp => real64

  ! We assume that kind=dp maps onto H5T_NATIVE_DOUBLE
  ! More generally, we could use h5kind_to_type(dp,H5_REAL_KIND)
  ! In the HDF5 Fortran 2003 API, a C pointer is used for the read/written dataset
  use iso_c_binding, only : c_ptr, c_loc
  use HDF5

  implicit none

  ! Public routines
  private
  public :: open_hdf5_file, close_hdf5_file
  public :: read_hdf5_dataset, write_hdf5_dataset
  public :: read_hdf5_attribute, write_hdf5_attribute

  ! File identifier and file_open flag. These are saved as module variables.
  integer(HID_T) :: file_id
  logical        :: file_open = .false.

  ! Generic procedures allowing a simple user interface

  interface write_hdf5_dataset
     module procedure write_hdf5_dataset_integer
     module procedure write_hdf5_dataset_real
  end interface write_hdf5_dataset

  interface read_hdf5_dataset
     module procedure read_hdf5_dataset_integer
     module procedure read_hdf5_dataset_real
  end interface read_hdf5_dataset

  interface write_hdf5_attribute
     module procedure write_hdf5_attribute_integer
     module procedure write_hdf5_attribute_real
     module procedure write_hdf5_attribute_string
  end interface write_hdf5_attribute

  interface read_hdf5_attribute
     module procedure read_hdf5_attribute_integer
     module procedure read_hdf5_attribute_real
     module procedure read_hdf5_attribute_string
  end interface read_hdf5_attribute

contains

  subroutine open_hdf5_file ( filename, readonly )
    implicit none
    character(*),      intent(in) :: filename ! Supplied file name
    logical, optional, intent(in) :: readonly ! Optionally indicates read only access

    integer        :: error     ! HDF5 error flag
    integer(HID_T) :: f_id      ! File id returned by HDF5 routine
    logical        :: read_only ! Local copy of readonly argument

    if ( file_open ) stop 'Attempt to open 2 HDF5 files at once'
    
    call h5open_f ( error ) ! Open Fortran interface

    read_only = .false.
    if ( present(readonly) ) read_only = readonly

    if ( read_only ) then
       write ( unit=output_unit, fmt='(a,t31,a)') 'Reading from HDF5 file',filename
       call h5fopen_f ( filename, H5F_ACC_RDONLY_F, f_id, error ) ! Open existing file, define file id
    else
       write ( unit=output_unit, fmt='(a,t31,a)') 'Writing to HDF5 file',filename
       call h5fcreate_f ( filename, H5F_ACC_TRUNC_F, f_id, error ) ! Create new file, define file id
    end if

    ! Store file_id and flag in module variables
    file_id   = f_id
    file_open = .true.

  end subroutine open_hdf5_file

  subroutine close_hdf5_file
    implicit none

    integer :: error ! HDF5 error flag

    write ( unit=output_unit, fmt='(a)') 'Closing HDF5 file'
    if ( .not. file_open ) stop 'Attempt to close an HDF5 file that is not open'
    
    call h5fclose_f ( file_id, error ) ! Close file
    call h5close_f  ( error ) ! Close Fortran interface

    ! Store flag in module variables
    file_open = .false.
    
  end subroutine close_hdf5_file

  subroutine write_hdf5_dataset_integer ( dataset_name, dataset )
    implicit none
    character(*),    intent(in) :: dataset_name ! Dataset name
    integer, target, intent(in) :: dataset(..)  ! Dataset values

    type(c_ptr)    :: dataset_ptr  ! C pointer to dataset
    integer(HID_T) :: dataset_id   ! Dataset identifier
    integer(HID_T) :: dataspace_id ! Dataspace identifier
    integer        :: error        ! Error flag

    write ( unit=output_unit, fmt='(a5,i1,a,t31,a)') 'Rank-', rank(dataset), ' dataset, integer:', dataset_name
    
    select rank ( dataset )
       rank(0)
          call h5screate_f ( H5S_SCALAR_F, dataspace_id, error )
       rank default
          call h5screate_simple_f ( rank(dataset), shape(dataset,kind=HSIZE_T), dataspace_id, error )
    end select

    call h5dcreate_f ( file_id, dataset_name, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error )
    dataset_ptr = c_loc(dataset)
    call h5dwrite_f  ( dataset_id, H5T_NATIVE_INTEGER, dataset_ptr, error )
    call h5dclose_f  ( dataset_id, error )
    call h5sclose_f  ( dataspace_id, error )

  end subroutine write_hdf5_dataset_integer

  subroutine read_hdf5_dataset_integer ( dataset_name, dataset )
    implicit none
    character(*),    intent(in)  :: dataset_name ! Dataset name
    integer, target, intent(out) :: dataset(..)  ! Dataset values

    type(c_ptr)    :: dataset_ptr ! C pointer to dataset
    integer(HID_T) :: dataset_id  ! Dataset identifier
    integer        :: error       ! Error flag

    write ( unit=output_unit, fmt='(a5,i1,a,t31,a)') 'Rank-', rank(dataset), ' dataset, integer:', dataset_name

    call h5dopen_f  ( file_id, dataset_name, dataset_id, error )
    dataset_ptr = c_loc(dataset)
    call h5dread_f  ( dataset_id, H5T_NATIVE_INTEGER, dataset_ptr, error )
    call h5dclose_f ( dataset_id, error )

  end subroutine read_hdf5_dataset_integer

  subroutine write_hdf5_dataset_real ( dataset_name, dataset )
    implicit none
    character(*),     intent(in) :: dataset_name ! Dataset name
    real(dp), target, intent(in) :: dataset(..)  ! Dataset values

    type(c_ptr)    :: dataset_ptr  ! C pointer to dataset
    integer(HID_T) :: dataset_id   ! Dataset identifier
    integer(HID_T) :: dataspace_id ! Dataspace identifier
    integer        :: error        ! Error flag

    write ( unit=output_unit, fmt='(a5,i1,a,t31,a)') 'Rank-', rank(dataset), ' dataset, real:', dataset_name

    select rank (dataset)
       rank(0)
          call h5screate_f ( H5S_SCALAR_F, dataspace_id, error )
       rank default
          call h5screate_simple_f ( rank(dataset), shape(dataset,kind=HSIZE_T), dataspace_id, error )
    end select
       
    call h5dcreate_f ( file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error )
    dataset_ptr = c_loc(dataset)
    call h5dwrite_f  ( dataset_id, H5T_NATIVE_DOUBLE, dataset_ptr, error )
    call h5dclose_f  ( dataset_id, error )
    call h5sclose_f  ( dataspace_id, error )

  end subroutine write_hdf5_dataset_real

  subroutine read_hdf5_dataset_real ( dataset_name, dataset )
    implicit none
    character(*),     intent(in)  :: dataset_name ! Dataset name
    real(dp), target, intent(out) :: dataset(..)  ! Dataset values

    type(c_ptr)    :: dataset_ptr ! C pointer to dataset
    integer(HID_T) :: dataset_id  ! Dataset identifier
    integer        :: error       ! Error flag

    write ( unit=output_unit, fmt='(a5,i1,a,t31,a)') 'Rank-', rank(dataset), ' dataset, real:', dataset_name

    call h5dopen_f  ( file_id, dataset_name, dataset_id, error )
    dataset_ptr = c_loc(dataset)
    call h5dread_f  ( dataset_id, H5T_NATIVE_DOUBLE, dataset_ptr, error )
    call h5dclose_f ( dataset_id, error )

  end subroutine read_hdf5_dataset_real

  subroutine write_hdf5_attribute_integer ( attribute_name, attribute )
    implicit none
    character(*),    intent(in) :: attribute_name ! Attribute name
    integer, target, intent(in) :: attribute(..)  ! Attribute values

    type(c_ptr)    :: attribute_ptr ! C pointer to attribute
    integer(HID_T) :: attribute_id  ! Attribute identifier
    integer(HID_T) :: attrspace_id  ! Attribute-space identifier
    integer        :: error         ! Error flag

    write ( unit=output_unit, fmt='(a5,i1,a,t31,a)') 'Rank-', rank(attribute), ' attribute, integer:', attribute_name

    select rank (attribute)
       rank(0)
          call h5screate_f ( H5S_SCALAR_F, attrspace_id, error )
       rank default
          call h5screate_simple_f ( rank(attribute), shape(attribute,kind=HSIZE_T), attrspace_id, error )
    end select

    call h5acreate_f ( file_id, attribute_name, H5T_NATIVE_INTEGER, attrspace_id, attribute_id, error )
    attribute_ptr = c_loc(attribute)
    call h5awrite_f ( attribute_id, H5T_NATIVE_INTEGER, attribute_ptr, error )
    call h5aclose_f ( attribute_id, error )
    call h5sclose_f ( attrspace_id, error )

  end subroutine write_hdf5_attribute_integer

  subroutine read_hdf5_attribute_integer ( attribute_name, attribute )
    implicit none
    character(*),    intent(in)  :: attribute_name ! Attribute name
    integer, target, intent(out) :: attribute(..)  ! Attribute value

    type(c_ptr)    :: attribute_ptr ! C pointer to attribute
    integer(HID_T) :: attribute_id  ! Attribute identifier
    integer        :: error         ! Error flag

    write ( unit=output_unit, fmt='(a5,i1,a,t31,a)') 'Rank-', rank(attribute), ' attribute, integer:', attribute_name

    call h5aopen_f  ( file_id, attribute_name, attribute_id, error )
    attribute_ptr = c_loc(attribute)
    call h5aread_f  ( attribute_id, H5T_NATIVE_INTEGER, attribute_ptr, error )
    call h5aclose_f ( attribute_id, error )

  end subroutine read_hdf5_attribute_integer

  subroutine write_hdf5_attribute_real ( attribute_name, attribute )
    implicit none
    character(*),     intent(in) :: attribute_name ! Attribute name
    real(dp), target, intent(in) :: attribute(..)  ! Attribute values

    type(c_ptr)    :: attribute_ptr ! C pointer to attribute
    integer(HID_T) :: attribute_id  ! Attribute identifier
    integer(HID_T) :: attrspace_id  ! Attribute-space identifier
    integer        :: error         ! Error flag

    write ( unit=output_unit, fmt='(a5,i1,a,t31,a)') 'Rank-', rank(attribute), ' attribute, real:', attribute_name

    select rank (attribute)
       rank(0)
          call h5screate_f ( H5S_SCALAR_F, attrspace_id, error )
       rank default
          call h5screate_simple_f ( rank(attribute), shape(attribute,kind=HSIZE_T), attrspace_id, error )
    end select

    call h5acreate_f ( file_id, attribute_name, H5T_NATIVE_DOUBLE, attrspace_id, attribute_id, error )
    attribute_ptr = c_loc(attribute)
    call h5awrite_f  ( attribute_id, H5T_NATIVE_DOUBLE, attribute_ptr, error )
    call h5aclose_f  ( attribute_id, error )
    call h5sclose_f  ( attrspace_id, error )

  end subroutine write_hdf5_attribute_real

  subroutine read_hdf5_attribute_real ( attribute_name, attribute )
    implicit none
    character(*),     intent(in)  :: attribute_name ! Attribute name
    real(dp), target, intent(out) :: attribute(..)  ! Attribute values

    type(c_ptr)    :: attribute_ptr ! C pointer to attribute
    integer(HID_T) :: attribute_id  ! Attribute identifier
    integer        :: error         ! Error flag

    write ( unit=output_unit, fmt='(a5,i1,a,t31,a)') 'Rank-', rank(attribute), ' attribute, real:', attribute_name

    call h5aopen_f  ( file_id, attribute_name, attribute_id, error )
    attribute_ptr = c_loc(attribute)
    call h5aread_f  ( attribute_id, H5T_NATIVE_DOUBLE, attribute_ptr, error )
    call h5aclose_f ( attribute_id, error )

  end subroutine read_hdf5_attribute_real

  subroutine write_hdf5_attribute_string ( attribute_name, attribute )
    implicit none
    character(*),         intent(in) :: attribute_name ! Attribute name
    character(*), target, intent(in) :: attribute      ! Attribute value

    type(c_ptr)    :: attribute_ptr ! C pointer to attribute
    integer(HID_T) :: attribute_id  ! Attribute identifier
    integer(HID_T) :: attrspace_id  ! Attribute-space identifier
    integer(HID_T) :: chartype_id   ! Character type identifier
    integer        :: error         ! Error flag

    write ( unit=output_unit, fmt='(a,t31,a)') 'Character string attribute:', attribute_name

    call h5tcopy_f ( H5T_NATIVE_CHARACTER, chartype_id, error )
    call h5tset_size_f ( chartype_id, len(attribute, kind=HSIZE_T), error )
    call h5screate_f ( H5S_SCALAR_F, attrspace_id, error )
    call h5acreate_f ( file_id, attribute_name, chartype_id, attrspace_id, attribute_id, error )
    attribute_ptr = c_loc(attribute)
    call h5awrite_f  ( attribute_id, chartype_id, attribute_ptr, error )
    call h5aclose_f  ( attribute_id, error )
    call h5sclose_f  ( attrspace_id, error )

  end subroutine write_hdf5_attribute_string

  subroutine read_hdf5_attribute_string ( attribute_name, attribute )
    implicit none
    character(*),         intent(in)  :: attribute_name ! Attribute name
    character(*), target, intent(out) :: attribute      ! Attribute value

    type(c_ptr)    :: attribute_ptr ! C pointer to attribute
    integer(HID_T) :: attribute_id  ! Attribute identifier
    integer(HID_T) :: chartype_id   ! Character type identifier
    integer        :: error         ! Error flag

    write ( unit=output_unit, fmt='(a,t31,a)') 'Character string attribute:', attribute_name

    call h5tcopy_f ( H5T_NATIVE_CHARACTER, chartype_id, error )
    call h5tset_size_f ( chartype_id, len(attribute, kind=HSIZE_T), error )
    call h5aopen_f ( file_id, attribute_name, attribute_id, error )
    attribute_ptr = c_loc(attribute)
    call h5aread_f  ( attribute_id, chartype_id, attribute_ptr, error )
    call h5aclose_f ( attribute_id, error )

  end subroutine read_hdf5_attribute_string

end module hdf5_module
