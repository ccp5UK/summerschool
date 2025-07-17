program block

!**********************************************************************
!     
!     blocking method for estimating standard error of mean
!     
!     (reference - Flyvbjerg and Petersen  JCP 91 (1989) 461)
!     
!     copyright daresbury laboratory
!     author - w. smith nov 1992
!     updated - j.purton 2017
!     
!**********************************************************************

    implicit none

    integer, parameter :: ndeg = 20, ndiv = 10000

    integer :: i, j, m, n, ncol, ndat
    real*8 :: ave, tim
    real*8 :: sample(3), d(ndiv), sigma(ndeg), gamma(ndeg)

    write(*,*)'CCP5 Summer School Statistical Error Program'

    write(*,*)"Select data for processing: 1=Eng,2=Tem,3=Prs"
    read(*,*) ncol

    ! open the statistical data file
      
    open(7,file="STA")
      
    ! read in statistical data (first column is time)

    do i = 1, ndiv

        read(7,*,end=100)tim, (sample(j),j=1,3)
        d(i) = sample(ncol)

    enddo

100 ndat=i-1

    close(7)
    write(*,*)"Number of data points read ",ndat

    ! calculate average

    ave = 0.d0

    do i = 1, ndat

        ave = ave + d(i)

    enddo

    ave = ave / dble(ndat)

    write(*,*)"Average value ", ave

    ! calculate errors using blocking method

    m = ndat

    do n = 1, ndeg

        sigma(n) = 0.d0

        if (n == 1) then

            do i = 1, m

                d(i) = d(i) - ave
                sigma(1) = sigma(1) + d(i)**2

            enddo

        else

            j = 0

            do i = 1, m

                j = j + 2
                d(i) = 0.5d0 * (d(j-1) + d(j))
                sigma(n) = sigma(n) + d(i)**2

            enddo

        endif

        sigma(n) = sqrt(sigma(n) / (dble(m) * dble(m-1)))
        if(m <= 3) go to 200
        m = m / 2

    enddo

200 continue
        
    ! calculate error bars of error estimates
      
    m = ndat
    do n = 1, ndeg

        gamma(n) = sigma(n) / sqrt(dble(2*(m-1)))
        if(m.le.3)go to 300
        m = m / 2

    enddo
300  continue
      
    write(*,*)'Estimated standard error and error bar'
    open(9,file="ERR")

    do i=1,min(n,ndeg)

        write(*,'(1p,e12.4,a,e11.4)')sigma(i)," +/-", gamma(i)
        write(9,'(i5,1p,2e12.4)')i, sigma(i), gamma(i)

    enddo

    close(9)

end




