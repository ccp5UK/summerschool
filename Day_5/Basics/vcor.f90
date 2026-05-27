program vcor
!*********************************************************************
!
!     dl_poly routine to calculate velocity autocorrelation function
!     for selected atoms
!
!     copyright daresbury laboratory 1995
!     author  w.smith jan 1995
!     sdk adaptation: january 1998
!
!*********************************************************************

    implicit none

    integer, parameter :: nvacf = 128

    integer :: i, j, l, m, lor, mor, iatms,natms
    integer :: nvm(nvacf), idv(nvacf)
    real*8 :: time, tzero, tstep, rnorm, vsum
    real*8, allocatable :: vacf(:,:), vcr0(:,:,:)
    real*8, allocatable :: vx(:), vy(:), vz(:)

    write(*,*)'CCP5 Summer School Velocity Autocorrelation Program'
    open (101,file="TRJ")
    read(101,*)natms, time
    close(101)
    allocate(vacf(natms,nvacf), vcr0(3,natms,nvacf))
    allocate(vx(natms), vy(natms), vz(natms) )

!   write control variables

    write(*,*)'Number of atoms in system      ', natms
    write(*,*)'Length of correlation arrays   ', nvacf

!   initialise velocity autocorrelation arrays

    lor = 0
    mor = 0
    do j = 1, nvacf

      nvm(j) = 0

      do i = 1, natms

        vacf(i,j) = 0.d0
        vcr0(1,i,j) = 0.d0
        vcr0(2,i,j) = 0.d0
        vcr0(3,i,j) = 0.d0

      enddo

  enddo

!   open trajectory file

open (7,file="TRJ")

!   process the trajectory file data

tzero = 0.d0
do while(.true.)

!     read trajectory file

  read(7,*,end=100)iatms, time
  tstep = time - tzero
  tzero = time

  if (natms /= iatms) then

    write(*,*)'error - incorrect number of atoms in TRJ file'
    stop

  endif

do i = 1, natms
  read(7,*) vx(i), vy(i), vz(i)
enddo

!     calculate velocity autocorrelation function

lor = min(lor+1,nvacf)
mor = mod(mor,nvacf) + 1
idv(mor) = 1

!     use new configuration as a time origin

do i = 1, natms

  vcr0(1,i,mor) = vx(i)
            vcr0(2,i,mor) = vy(i)
            vcr0(3,i,mor) = vz(i)

        enddo

!     calculate velocity autocorrelation
        do l = 1, lor

            m = idv(l)
            idv(l) = m + 1
            nvm(m) = nvm(m) + 1

            do i = 1, natms

                vacf(i,m) = vacf(i,m) + vx(i) * vcr0(1,i,l) + vy(i) * vcr0(2,i,l) + vz(i) * vcr0(3,i,l)

            enddo

        enddo

    enddo

100 continue

!     normalise velocity autocorrelation function

    do i = 1, natms

        rnorm = dble(nvm(1)) / vacf(i,1)

        do j = 1, lor

          vacf(i,j) = rnorm * vacf(i,j) / dble(nvm(j))

        enddo

    enddo

!     print out final velocity correlation function

    open(8,file="VAF")
    do j = 1, lor

        vsum = 0.0d0

        do i = 1, natms

          vsum = vsum + vacf(i,j)

        enddo

        vsum = vsum / dble(natms)
        time = tstep * dble(j)
        write(8,'(2f16.7)')time, vsum

    enddo

    close(7)
    close(9)
    write(*,*)'Job done'
end
