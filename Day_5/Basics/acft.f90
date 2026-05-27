program acft

!***********************************************************************
!     
!     dl_poly program for calculating correlation functions
!     using the fast fourier transform strategy
!     
!     copyright - daresbury laboratory 1992
!     author    - w. smith may 1992.
!     
!     note - this program includes a naff fft routine purely for
!     portability and verification reasons. users are strongly
!     recommended to substitute a machine specific fft, which
!     should vastly improve performance.
!     
!     note - you will probably need to change some of the array
!     dimensions using the following parameters. 
!     
!     note - ndiv must be a power of 2. ndiv3 will change according to
!     the fft routine you use. leave ndiv2 as it is.
!     
!***********************************************************************

    implicit none

    integer, parameter :: ndiv = 16384, ndiv2 = 2 * ndiv, ndiv3 = 2 * ndiv
      
    integer :: i, j, ncol, ndat
    integer :: key(ndiv2)
    real*8 :: dtim, rnrm, time, ave
    real*8 :: sample(3), u(ndiv), t(ndiv)
    complex*16 :: g(ndiv2), c(ndiv2), w(ndiv3)
      
    write(*,*)'CCP5 Summer School Velocity Autocorrelation Program'
      
    write(*,*)'Select a data number 1-3: (Egy, Tem, Prs)'
    read(*,*)ncol

    ! zero the complex data array
         
    do i=1,ndiv2
        g(i) = (0.d0,0.d0)
    enddo
         
    ! open the statistical data file
         
    open(7,file="STA")
         
    ! now read the required column (note first column is time)
         
    do i=1,ndiv

        read(7,*,end=100)t(i), (sample(j),j=1,3)
        g(i) = sample(ncol)

    enddo
100 ndat = i - 1

    close(7)

    ! time interval
         
    dtim = t(2) - t(1)
         
    ! calculate average
         
    ave = 0.d0
    do i = 1, ndat
        ave = ave + real(g(i))
    enddo
    ave = ave / dble(ndat)
         
    ! subtract average
         
    do i = 1, ndat
        g(i) = g(i) - ave
    enddo
    write(*,*)'Average :', ave
      
    ! initialise the fft work arrays
      
    call fft(1, 1, ndiv2, key, g, w, g)
      
    ! fourier transform the data
      
    call fft(0, -1, ndiv2, key, g, w, g)
      
    ! calculate correlation as product in frequency domain
      
    do i = 1, ndiv2

        c(i) = g(i) * conjg(g(i))

    enddo
      
    ! inverse fourier transform
      
    call fft(0, 1, ndiv2, key, c, w, c)
      
    ! normalise the correlation function
      
    rnrm = 1.d0 / dble(ndiv2)
    u(1) = rnrm * real(c(1)) / dble(ndat)

    write(*,*)'Normalisation :',u(1)

    do i = 2, ndat
        u(i) = (rnrm * real(c(i)) / dble(ndat+1-i)) / u(1)
    enddo

    u(1) = 1.d0

    ! open correlation function file
      
    open(8,file="COR")

    ! write out correlation function
      
    do i=1,ndat

        time = dble(i-1) * dtim
        write(8,'(f10.6,1p,e15.7)')time, u(i)

    enddo

    close(8)
    write(*,*)'Job done'

end
      
subroutine fft(ind, isw, ndiv, key, aaa, wfft, bbb)
!***********************************************************************
!     
!     fast fourier transform routine
!     
!***********************************************************************
      
    implicit none
      
    integer :: i, j, l, nt, nu, ind, isw, iii, jjj, kkk, kk1, k12, jj2
    integer :: ndiv, np1, np2, nu1
    integer :: key(ndiv)
    real*8 :: tpi, tpn, arg
    complex*16 :: aaa(ndiv), bbb(ndiv), wfft(ndiv), ttt, cmplx

    logical :: check
    data tpi/6.2831853072d0/
      
    ! check that array is of suitable length

    nt = 1
    check = .true.

    do i = 1, 20

        nt = 2 * nt

        if (nt == ndiv) then

            check = .false.
            nu = i

        endif

    enddo

    if (check) then

         write(*,*)'error - number of points not a power of two'
         stop

    endif
      
    if (ind > 0) then

    ! set reverse bit address array

        do kkk=1,ndiv

            iii = 0
            jjj = kkk - 1

            do j = 1, nu

               jj2 = jjj / 2
               iii = 2 * (iii - jj2) + jjj
               jjj = jj2

            enddo

            key(kkk) = iii + 1

        enddo

        ! initialise complex exponential factors
         
        tpn = tpi / dble(ndiv)
        arg = 0.d0
        np1 = ndiv + 1
        np2 = ndiv / 2
        wfft(1) = (1.d0,0.d0)

        do i = 1, np2

            arg = tpn * dble(i)
            wfft(i+1) = cmplx(cos(arg),sin(arg))
            wfft(np1-i) = conjg(wfft(i+1))

        enddo

        return

    endif

    ! take conjugate of exponentials if required
      
    if (isw < 0) then

        do i = 1, ndiv

             wfft(i) = conjg(wfft(i))

        enddo

    endif
      
    ! take copy input array
      
    do i = 1, ndiv

        bbb(i)=aaa(i)

    enddo

    ! perform fourier transform
      
    kkk = 0
    nu1 = nu - 1
    np2 = ndiv / 2

    do l = 1, nu

100     do i = 1, np2

            iii = key(kkk/2**nu1+1)
            kk1 = kkk + 1
            k12 = kk1 + np2
            ttt = bbb(k12) * wfft(iii)
            bbb(k12) = bbb(kk1) - ttt
            bbb(kk1) = bbb(kk1) + ttt
            kkk = kkk + 1

        enddo

        kkk = kkk + np2

        if (kkk < ndiv) go to 100

        kkk = 0
        nu1 = nu1 - 1
        np2 = np2 / 2

    enddo

    ! unscramble the fft using bit address array
      
    do kkk = 1, ndiv

        iii = key(kkk)

        if (iii > kkk) then

            ttt = bbb(kkk)
            bbb(kkk) = bbb(iii)
            bbb(iii) = ttt

        endif

    enddo

    ! restore exponentials to unconjugated values if necessary
      
    if (isw < 0) then

        do i = 1, ndiv

            wfft(i) = conjg(wfft(i))

        enddo

    endif
      
    return

end
