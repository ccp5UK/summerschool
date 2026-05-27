program auto
!***********************************************************************
!     
!     dl_poly program for calculating autocorrelation functions
!     
!     copyright - daresbury laboratory 1992
!     author    - w. smith may 1992.
!     
!***********************************************************************
      
    implicit none

    integer, parameter :: ndiv = 10000
      
    integer :: i, j, ncol, ndat
    real*8 :: dtim, time, sum
    real*8 :: c(ndiv), u(ndiv), t(ndiv), sample(3)

    write(*,*)'CCP5 Summer School Autocorrelation Program'

    write(*,*)'Select a data number 1-3: (Egy, Tem, Prs)'
    read(*,*) ncol

    ! open the statistical data file

    open(7,file="STA")

    ! now read the required column (note first column is time)
      
    do i = 1, ndiv

        read(7,*,end=100) t(i), (sample(j),j=1,3)
        u(i) = sample(ncol)

    enddo
      
100 ndat = i - 1
    close(7)

    ! time interval
      
    dtim = t(2) - t(1)
      
    ! calculate average
      
    sum = 0.d0
    do i = 1, ndat

        sum = sum + u(i)

    enddo

    sum = sum / dble(ndat)
    write(*,*)'Average :', sum

    ! subtract average
      
    do i = 1, ndat

        u(i) = u(i) - sum

    enddo

    ! calculate correlation function
      
    do i = 1, ndat

        c(i) = 0.d0

        do j = 1, ndat + 1 - i

          c(i) = c(i) + u(j) * u(i+j-1)

        enddo

    enddo
      
    ! normalise the correlation function

    u(1) = c(1) / dble(ndat)

    write(*,*)'Normalisation :', u(1)
    do i = 2, ndat

        u(i) = (c(i) / dble(ndat+1-i)) / u(1)

    enddo

    u(1) = 1.d0

    ! open correlation function file
      
    open(8,file="COR")

    ! write out correlation function
      
    do i = 1, ndat

        time = dble(i-1) * dtim
        write(8,'(f10.6,1p,e15.7)')time, u(i)

    enddo

    close(8)
    write(*,*)'Job done'

end
