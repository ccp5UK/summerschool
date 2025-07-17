!    Copyright (c) 2014-2019 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module mLife
  implicit none

  private
  public :: initializeBoard
  public :: printBoard
  public :: evolveBoard
  public :: census
  public :: swapHalos

contains

  subroutine initializeBoard(a,n,m)
    integer, intent(in)    :: n,m
    integer, intent(inout),allocatable :: a(:,:)

    integer      :: i,j
    real(kind=8) :: x
    a(0,:)=0
    a(n+1,:)=0
    a(:,0)=0
    a(:,m+1)=0
    do j=1,m
      do i=1,n
        call random_number(x)
        if (x<0.5_8) then
          a(i,j)=1
        else
          a(i,j)=0
        endif
      enddo
    enddo
  end subroutine initializeBoard

  subroutine printBoard(ffile,a,n,m)
    integer, intent(in)    :: n,m
    integer, intent(in),allocatable :: a(:,:)
    character(len=*), intent(in)  :: ffile

    integer :: i
    integer :: u

    open(newunit=u,file=trim(ffile), action="write",status="unknown")

    do i=1,n
        write(u,'(*(i0,1x))')a(i,:)
      write(u,*)
    enddo
    close(u)
  end subroutine printBoard

  subroutine evolveBoard(a,n,m)
    integer, intent(in)    :: n,m
    integer, intent(inout),allocatable :: a(:,:)

    integer :: i,j,ls,b(n,m)

    do j=1,m
      do i=1,n
        ls=sum(a(i-1:i+1,j-1))+sum(a(i-1:i+1,j+1))+a(i-1,j)+a(i+1,j)
        if ((ls < 2) .or. (ls > 3)) b(i,j)=0
        if (ls==3) b(i,j)=1
        if (ls==2) b(i,j)=a(i,j)
      enddo
    enddo
    a(1:n,1:m)=b(1:n,1:m)
  end subroutine evolveBoard

  integer function census(a,n,m)
    integer, intent(in)    :: n,m
    integer, intent(in),allocatable :: a(:,:)
    census=sum(a(1:n,1:m))
  end function census

  subroutine swapHalos(a,n,m)
    integer, intent(in)    :: n,m
    integer, intent(inout),allocatable :: a(:,:)

    ! west to east
    a(:,0)=a(:,m)
    a(:,m+1)=a(:,1)
    a(n+1,:)=a(1,:)
    a(0,:)=a(n,:)
  end subroutine swapHalos

end module mLife

program life
  use mLife
  use iso_fortran_env , only : error_unit
  implicit none
  integer, allocatable        :: a(:,:)

  integer           :: n,m
  integer           :: ierror
  character(len=32) :: arg
  integer           :: i,s,lpopulation, nGenerations
  if (command_argument_count() /=3) then
    write(*,'(a)')"Wrong usage!!!"
    write(*,'(a)')"Correct usage: "
    call get_command_argument(0, arg)
    write(*,'(a)')trim(arg)//" <n> <m> <generations>"
    write(*,'(a)')trim(arg)//" 10 12 10"
    error stop -1
  endif

  call get_command_argument(1, arg)
  read(arg,*) n
  call get_command_argument(2, arg)
  read(arg,*) m
  call get_command_argument(3, arg)
  read(arg,*) nGenerations

  write(*,'(a)',advance="no")"Running: "
  do i=0,command_argument_count()
    call get_command_argument(i, arg)
    write(*,'(a,1x)',advance="no") trim(arg)
  enddo
  write(*,*)
! allocate the memory... distributed per image
  allocate(a(0:n+1,0:m+1),stat=ierror)
  if (ierror/=0) then
    write(error_unit,'(a)')"Error -2: Failed to allocate memory"
    error stop
  end if
! initialize the board
  call initializeBoard(a,n,m)
  call swapHalos(a,n,m)
  print *, "Initial population"
  call printBoard("life_init.dat",a,n,m)
  lpopulation=census(a,n,m)
  open(101,file="life.dat",status="unknown",action="write")
  write(101,'(a)')"#Generation | Population "
  write(101,'(i0,1x,i0)')0,lpopulation
  do i=1,nGenerations
    call evolveBoard(a,n,m)
    lpopulation=census(a,n,m)
    write(101,'(i0,1x,i0)')i,lpopulation
    call swapHalos(a,n,m)
  end do
  close(101)
  print *, "Final population"
  call printBoard("life_final.dat",a,n,m)
  deallocate(a)
end program life

