!---------------------------------------------------------------------
!     
!     program to compare leapfrog algorithm with alternative
!     integration algorithms:
!     
!     - original Verlet algorithm
!     - velocity Verlet algorithm
!     - Rahman algorithm
!     
!     copyright daresbury laboratory
!     author w.smith 2003
!     updated j.purton 2017
!     
!---------------------------------------------------------------------

      implicit none

      integer :: i, j, icyc, nsteps, nprnt, numacc, key
      real*8 tstep, pot, frc, x, y, dev, eng, ena, time
      real*8 stpval, sclnv1, sclnv2, avgeng, av2eng, avgena, av2ena
      real*8 xxx, yyy, vxx, vyy, fxx, fyy
      real*8 xxa, yya, xxb, yyb, xxc, yyc, vxa, vya, vxc, vyc, fxa, fya, fxc, fyc

      !     define potential and force

      pot(x,y) = 0.5d0 * x * x * y * y * (x * x + y * y)
      frc(x,y) = -x * y * (2.d0 * x * x * y + y * y * y)

      !     open output file

      open(7,file="OUT")

      icyc=2
      !     note: icyc must be >= 1
      numacc = 0
      read(*,*) key
      read(*,*) nsteps
      read(*,*) nprnt
      read(*,*) tstep

      avgeng = 0.d0
      avgena = 0.d0

      !     velocity verlet initial position, velocity and force

      xxx = 1.d0
      yyy = 1.d0
      vxx = 1.d0
      vyy = 1.d0 / 3.d0

      fxx = frc(xxx,yyy)
      fyy = frc(yyy,xxx)

      if (key == 1) then
        
      !     leapfrog initial position and velocity

          xxa = xxx
          yya = yyy
          vxa = vxx - 0.5d0 * tstep * fxx
          vya = vyy - 0.5d0 * tstep * fyy
        
      else if (key == 2) then
        
      !     verlet initial positions

          xxa = xxx
          yya = yyy
          vxa = vxx - 0.5d0 * tstep * fxx
          vya = vyy - 0.5d0 * tstep * fyy
          xxb = xxx - tstep * vxa
          yyb = yyy - tstep * vya

      else if (key >= 3) then

      !     rahman initial position, velocity and force 

          xxa = xxx
          yya = yyy
          vxa = vxx
          vya = vyy
          fxa = fxx
          fya = fyy
          vxc = vxx - 0.5d0 * tstep * fxx
          vyc = vyy - 0.5d0 * tstep * fyy
          xxb = xxx - tstep * vxc
          yyb = yyy - tstep * vyc
        
      endif

      do i = 0, nsteps

      !     velocity verlet algorithm

          vxx = vxx + 0.5d0 * tstep * fxx
          vyy = vyy + 0.5d0 * tstep * fyy
          xxx = xxx + tstep * vxx
          yyy = yyy + tstep * vyy
          fxx = frc(xxx,yyy)
          fyy = frc(yyy,xxx)
          vxx = vxx + 0.5d0 * tstep * fxx
          vyy = vyy + 0.5d0 * tstep * fyy

          if (key == 1) then

      !     leapfrog algorithm

              fxa = frc(xxa,yya)
              fya = frc(yya,xxa)
              vxa = vxa + tstep * fxa
              vya = vya + tstep * fya
              xxa = xxa + tstep * vxa
              yya = yya + tstep * vya

          else if (key == 2) then

      !     original verlet algorithm

              fxa = frc(xxa,yya)
              fya = frc(yya,xxa)
              xxc = xxa
              yyc = yya
              xxa = 2.d0 * xxa - xxb + fxa * tstep**2
              yya = 2.d0 * yya - yyb + fya * tstep**2
              xxb = xxc
              yyb = yyc

          else if (key >= 3) then
          
      !     rahman algorithm predictor
              xxc = xxb + 2.d0 * tstep * vxa
              yyc = yyb + 2.d0 * tstep * vya

      !     rahman algorithm corrector
              do j = 1, icyc

                  fxc = frc(xxc,yyc)
                  fyc = frc(yyc,xxc)
                  vxc = vxa + 0.5d0 * tstep * (fxc + fxa)
                  vyc = vya + 0.5d0 * tstep * (fyc + fya)
                  xxc = xxa + 0.5d0 * tstep * (vxc + vxa)
                  yyc = yya + 0.5d0 * tstep * (vyc + vya)

              enddo

       !     rahman final position
              xxb = xxa
              yyb = yya
              xxa = xxc
              yya = yyc
              vxa = vxc
              vya = vyc
              fxa = fxc
              fya = fyc

          endif

       !     calculate deviation

          if (mod(i,nprnt) == 0) then

       !     calculate energy and related statistics

              time = dble(i) * tstep
              if (key == 4) then

                  numacc = numacc + 1
                  sclnv2 = 1.d0 / dble(numacc)
                  sclnv1 = dble(numacc - 1) / dble(numacc)
                  eng = 0.5d0 * (vxx * vxx + vyy * vyy) + pot(xxx,yyy)
                  ena = 0.5d0 * (vxa * vxa + vya * vya) + pot(xxa,yya)
                  write(7,'(8e18.10)') time, eng, ena
                  av2eng = sclnv1 * (av2eng + sclnv2 * (eng - avgeng)**2)
                  av2ena = sclnv1 * (av2ena + sclnv2 * (ena - avgena)**2)
                  avgeng = sclnv1 * avgeng + sclnv2 * eng
                  avgena = sclnv1 * avgena + sclnv2 * ena
          
              else

                  dev = sqrt((xxx - xxa)**2 + (yyy - yya)**2)
                  write(7,'(8e18.10)')time, xxx, yyy, xxa, yya, dev

              endif

          endif

       !     time reversal switch

          if (key == 5) then

              if(i == nsteps/2) then
            
                  vxx = -vxx
                  vyy = -vyy
                  vxa = -vxa
                  vya = -vya
            
              endif

          endif

       !     end of dynamics loop

      enddo

       !     write energy statistics

      if(key.eq.4)then

          write(*,*)'v. verlet energy and fluctuation'
          write(*,'(2e18.10)')avgeng, sqrt(av2eng)

          write(*,*)'rahman energy and fluctuation'
          write(*,'(2e18.10)')avgena, sqrt(av2ena)

      endif

       !     close output file

      close(7)

      end
