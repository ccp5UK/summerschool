program msdp
!*********************************************************************
!
!     dl_poly program to calculate mean square displacement
!
!     copyright daresbury laboratory 1996
!     author  w.smith jan 1996
!     sdk adaptation: january 1998
!
!*********************************************************************

    implicit none

    integer, parameter :: nmsd = 512

    integer :: i, j, m, lsr, msr, iatms,natms
    integer :: msm(nmsd), imd(nmsd)

    real*8 :: time, tzero, tstep, dx, dy, dz, box, rmsx, rmsy, rmsz, rr2
    real*8,allocatable :: x(:), y(:), z(:)
    real*8,allocatable :: x0(:), y0(:), z0(:)
    real*8,allocatable :: ax(:), ay(:), az(:)
    real*8,allocatable :: msd(:,:), msd0(:,:,:)

    logical :: switch

    write(*,*)'CCP5 Summer School Mean Squared Displacement Program'

    open (101,file="TRJ")
    read(101,*)natms, time, box
    close(101)
    allocate( x(natms), y(natms), z(natms) )
    allocate(x0(natms),y0(natms), z0(natms), ax(natms), ay(natms),az(natms))
    allocate( msd(natms,nmsd),msd0(3,natms,nmsd))

    switch = .false.

    ! write control variables

    write(*,*)'Number of atoms in system      ', natms
    write(*,*)'Length of correlation arrays   ', nmsd

    ! initialise msd arrays

    lsr = 0
    msr = 0
    do i = 1, natms

      ax(i) = 0.d0
      ay(i) = 0.d0
      az(i) = 0.d0

    enddo

  do j=1,nmsd

    msm(j) = 0

    do i = 1, natms

      msd(i,j) = 0.d0
      msd0(1,i,j) = 0.d0
      msd0(2,i,j) = 0.d0
      msd0(3,i,j) = 0.d0

    enddo

enddo

 ! open trajectory file

open (7,file="TRJ")

 ! process the trajectory file data

tzero=0.d0
do while(.true.)

  ! read trajectory file

  read(7,*,end=100)iatms, time, box
  tstep = time - tzero
  tzero = time

  if(natms /= iatms)then

            write(*,*)'error - incorrect number of atoms in TRJ file'
            stop

        endif

        do i=1, natms

            read(7,*) x(i), y(i), z(i)

        enddo

    ! accumulate incremental distances

        if (switch) then

            do i = 1, natms

                dx = x(i) - x0(i)
                dy = y(i) - y0(i)
                dz = z(i) - z0(i)
                dx = dx - box * nint(dx / box)
                dy = dy - box * nint(dy / box)
                dz = dz - box * nint(dz / box)
                ax(i) = ax(i) + dx
                ay(i) = ay(i) + dy
                az(i) = az(i) + dz

            enddo

        endif

        do i = 1, natms

            x0(i) = x(i)
            y0(i) = y(i)
            z0(i) = z(i)

        enddo

    ! calculate mean square displacement

        if (switch) then

            lsr = min(lsr+1,nmsd)
            msr = mod(msr,nmsd) + 1
            imd(msr) = 1

            do i = 1, natms

                msd0(1,i,msr) = 0.d0
                msd0(2,i,msr) = 0.d0
                msd0(3,i,msr) = 0.d0

            enddo

        endif

        do j = 1, lsr

            m = imd(j)
            imd(j) = m + 1
            msm(m) = msm(m) + 1

            do i = 1, natms

                rmsx = msd0(1,i,j) + ax(i)
                rmsy = msd0(2,i,j) + ay(i)
                rmsz = msd0(3,i,j) + az(i)
                msd(i,m) = msd(i,m) + rmsx**2 + rmsy**2 + rmsz**2
                msd0(1,i,j) = rmsx
                msd0(2,i,j) = rmsy
                msd0(3,i,j) = rmsz

            enddo

        enddo

        do i = 1, natms

            ax(i) = 0.d0
            ay(i) = 0.d0
            az(i) = 0.d0

        enddo

        switch=.true.

    enddo

100 continue

    ! normalise and print mean square displacement

    open(8,file="MSD")

    do j = 1, lsr

        rr2 = 0.d0

        do i = 1, natms

          rr2 = rr2 + msd(i,j)

        enddo

        rr2 = rr2 / (dble(natms) * dble(msm(j)))
        time = tstep * dble(j)
        write(8,'(2f16.7)')time, rr2

    enddo

    close(7)
    close(8)
    write(*,*)'Job done'

end
