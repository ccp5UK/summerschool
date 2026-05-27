program sfac

!*********************************************************************
!
!     dl_poly/sdk routine to read a DL_POLY RDFDAT file
!
!     copyright - daresbury laboratory
!     author    - w.smith january 1998
!
!*********************************************************************

    implicit none

    integer, parameter :: ngrid=128

    integer :: i,j,isw,npts

    real*8 :: delr
    real*8 :: rrr(ngrid), gor(ngrid)

    write(*,*)'CCP5 Summer School Structure Factor Program'

    ! open RDF data file

    open(7,file="RDF")

    do i = 1, ngrid

        gor(i) = 0.d0

    enddo

    do i = 1, ngrid

        read(7,*,end=100)rrr(i), gor(i)
        gor(i) = gor(i) - 1.d0

    enddo
100 continue

    npts = i - 1
    close(7)

    do i = 1, npts - 1

        rrr(i) = 0.5 * (rrr(i) + rrr(i+1))
        gor(i) = 0.5 * (gor(i) + gor(i+1))

    enddo

    npts = npts - 1
    delr = (rrr(npts) - rrr(1)) / dble(npts-1)

    ! calculate radial Fourier transform of RDF

    isw = 1
    call radfft(isw, ngrid, delr, rrr, gor)

    ! print out SOK plot file

    open(8,file="SOK")

    do i = 1, ngrid

        write(8,'(1x,1p,2e15.7)')rrr(i), gor(i)

    enddo

    close(8)
    write(*,*)"Job done"

end
subroutine radfft(isw, nnn, delr, aaa, bbb)
!***********************************************************************
!
!     dl_poly 3D radial fourier transform routine using lado's method
!     reference: j. comput. phys. 8 (1971) 417
!
!     copyright daresbury laboratory 1994
!
!     author w smith
!
!     note: first data point is i*delr not 0
!
!***********************************************************************

    implicit none

    real*8, parameter :: pi=3.141592653589793d0

    integer :: isw, nnn, i, j
    real*8 :: delr, sw, delk
    real*8 :: aaa(nnn), bbb(nnn)

    ! perform fourier transform

    sw = pi * dble(isw) / dble(nnn)

    do j = 1, nnn

        aaa(j) = 0.d0

        do i = 1, nnn

            aaa(j) = dble(i) * bbb(i) * sin(sw * dble(j) * dble(i)) + aaa(j)

        enddo

        aaa(j) = (4.d0 * dble(nnn) * delr**3 / dble(j)) * aaa(j)

    enddo

    delk = pi / (delr * dble(nnn))
    do i = 1, nnn

        bbb(i) = aaa(i)
        aaa(i) = delk * dble(i)

    enddo

    return
end
