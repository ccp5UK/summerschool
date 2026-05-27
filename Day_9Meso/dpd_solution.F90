MODULE variables
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND (15, 307)
      REAL(KIND=dp), ALLOCATABLE ::   x(:),  y(:),  z(:)
      REAL(KIND=dp), ALLOCATABLE ::  vx(:), vy(:), vz(:)
      REAL(KIND=dp), ALLOCATABLE ::  fx(:), fy(:), fz(:)
#ifdef FAST
      REAL(KIND=dp) :: wdthlnk
      INTEGER, ALLOCATABLE :: lct(:), link(:)
      INTEGER :: nlink
#endif
      REAL(KIND=dp) :: boxside, poteng, virial, aij, bij, gamma, sigma, fac12
      INTEGER :: nbeads
END MODULE variables

PROGRAM dpd

!**********************************************************************
!     
!     CCP5 Summer School Dissipative Particle Dynamics Program
!     
!     copyright:  STFC Daresbury Laboratory
!     author:     M. A. Seaton, 2021
!     
!     program:    single component systems (optionally using link-cell
!                 method)
!
!***********************************************************************

      USE variables
      IMPLICIT none

      LOGICAL :: ltraj, lpos
      INTEGER :: i, j, k, l, n, istep, nstep, nequil, nav, isample, mbeads, nrep
      REAL(KIND=dp) :: displ, mx, my, mz, scl, tstephalf, fac, time, timesim
      REAL(KIND=dp) :: toteng, totengav, kineng, temp, tempav, tempsys, tstep, volm
      REAL(KIND=dp) :: press, pressav 

      WRITE (*,*)'CCP5 Summer School Dissipative Particle Dynamics Program'

!     read control data from input file

      READ(*,*)nstep
      READ(*,*)nequil
      READ(*,*)isample
      READ(*,*)nbeads
      READ(*,*)volm
      READ(*,*)tempsys
      READ(*,*)aij
      READ(*,*)bij
      READ(*,*)gamma
      READ(*,*)tstep
      READ(*,*)ltraj
      READ(*,*)lpos

      WRITE(*,*)'Number of time steps          ',nstep
      WRITE(*,*)'Number of equilibration steps ',nequil
      WRITE(*,*)'Properties sampling interval  ',isample
      WRITE(*,*)'Number of beads               ',nbeads
      WRITE(*,*)'System volume                 ',volm
      WRITE(*,*)'System temperature            ',tempsys
      WRITE(*,*)'Energy parameter              ',aij
      WRITE(*,*)'Particle diameter             ',bij
      WRITE(*,*)'Drag coefficient              ',gamma
      WRITE(*,*)'Simulation time step          ',tstep
      WRITE(*,*)'Trajectory file switch        ',ltraj
      WRITE(*,*)'Trajectory file positions     ',lpos

!     set up simulation

      nav = 0
      toteng = 0.0_dp
      totengav = 0.0_dp
      pressav = 0.0_dp
      tempav = 0.0_dp
      fac = 1.5_dp
      fac12 = 2.0_dp * SQRT(3.0_dp)
      tstephalf = tstep * 0.5_dp
      boxside = volm**(1.0_dp/3.0_dp)
      sigma = SQRT (2.0_dp * gamma * tempsys / tstep)
      nrep = MAX(1,nstep/20)

#ifdef FAST
!     set up link-cell structure (rcut = 1.0_dp)

      nlink = INT (boxside)
      IF (nlink==0) nlink = 1
      wdthlnk = boxside / DBLE (nlink)
#endif

!     allocate arrays for beads

      ALLOCATE (x (nbeads), y (nbeads), z (nbeads))
      ALLOCATE (vx (nbeads), vy (nbeads), vz (nbeads))
      ALLOCATE (fx (nbeads), fy (nbeads), fz (nbeads))
#ifdef FAST
      ALLOCATE (lct (nbeads), link (nbeads))
#endif

!     set initial positions of beads (in cubic lattice)

      mbeads = nbeads**(1.0_dp/3.0_dp)
      DO WHILE (mbeads*mbeads*mbeads<nbeads)
        mbeads = mbeads + 1
      END DO
      displ = boxside / REAL (mbeads+1, KIND=dp)

      DO k = 0, mbeads - 1
        DO j = 0, mbeads - 1 
          DO i = 0, mbeads - 1
            n = 1 + i + mbeads * (j + k * mbeads)
            IF (n<=nbeads) THEN
              x(n) = (REAL(i, KIND=dp) + 0.5_dp) * displ - 0.5_dp * boxside
              y(n) = (REAL(j, KIND=dp) + 0.5_dp) * displ - 0.5_dp * boxside
              z(n) = (REAL(k, KIND=dp) + 0.5_dp) * displ - 0.5_dp * boxside
            END IF
          END DO
        END DO
      END DO

!     set initial velocities

      mx = 0.0_dp
      my = 0.0_dp
      mz = 0.0_dp
      CALL RANDOM_NUMBER(vx)
      CALL RANDOM_NUMBER(vy)
      CALL RANDOM_NUMBER(vz)
      DO k = 1, nbeads
        vx(k) = vx(k) - 0.5_dp
        vy(k) = vy(k) - 0.5_dp
        vz(k) = vz(k) - 0.5_dp
        mx = mx + vx(k)
        my = my + vy(k)
        mz = mz + vz(k)
      END DO
      mx = mx / REAL(nbeads, KIND=dp)
      my = my / REAL(nbeads, KIND=dp)
      mz = mz / REAL(nbeads, KIND=dp)

!     remove net system momentum

      kineng = 0.0_dp
      DO k = 1, nbeads
        vx(k) = vx(k) - mx
        vy(k) = vy(k) - my
        vz(k) = vz(k) - mz
        kineng = kineng + vx(k) * vx(k) + vy(k) * vy(k) + vz(k) * vz(k)
      END DO
      kineng = 0.5_dp * kineng

!     temperature scaling

      temp = kineng / (fac * REAL(nbeads, KIND=dp))
      scl = SQRT (temp/tempsys)
      vx = scl * vx
      vy = scl * vy
      vz = scl * vz

!     open statistics file

      OPEN (9, FILE = "STATS")

!     open trajectory file

      IF (ltraj) OPEN (11, FILE = "TRAJECT")

!     calculate initial forces

#ifdef FAST
      CALL linkcell
#endif
      CALL dpd_forces

!     initialise clock for timing simulation

      CALL timecheck (timesim)

!     velocity verlet algorithm

      DO istep = 1, nstep

        poteng = 0.0_dp
        kineng = 0.0_dp
        virial = 0.0_dp
        
!     first half step motion
        
        DO k = 1, nbeads
          vx(k) = vx(k) + tstephalf * fx(k)
          vy(k) = vy(k) + tstephalf * fy(k)
          vz(k) = vz(k) + tstephalf * fz(k)
          x(k) = x(k) + tstep * vx(k)
          y(k) = y(k) + tstep * vy(k)
          z(k) = z(k) + tstep * vz(k)
        END DO
        
!     restore periodic boundary
        
        DO k = 1, nbeads
          x(k) = x(k) - boxside * NINT (x(k)/boxside)
          y(k) = y(k) - boxside * NINT (y(k)/boxside)
          z(k) = z(k) - boxside * NINT (z(k)/boxside)
        END DO

#ifdef FAST
!     construct link-cell lists

        CALL linkcell

#endif
!     calculate new forces

        CALL dpd_forces

!     second half step motion
        
        DO k = 1, nbeads
          vx(k) = vx(k) + tstephalf * fx(k)
          vy(k) = vy(k) + tstephalf * fy(k)
          vz(k) = vz(k) + tstephalf * fz(k)
          kineng = kineng + vx(k) * vx(k) + vy(k) * vy(k) + vz(k) * vz(k)
        END DO
        kineng = 0.5_dp * kineng

!     system properties

        temp = kineng / (fac * nbeads)
        toteng = kineng + poteng
        press = (2.0_dp * kineng - virial) / (3.0_dp * volm)
        IF (istep>nequil) THEN
          nav = nav + 1
          tempav = ((tempav * (nav-1)) / nav) + (temp/nav)
          totengav = ((totengav * (nav-1)) / nav) + (toteng/nav)
          pressav = ((pressav * (nav-1)) / nav) + (press/nav)
        END IF
        
!     checkpoint report

        time = tstep * REAL(istep, KIND=dp)
        IF (MOD (istep, nrep)==0) WRITE (*,'(1x,1p,4e16.8)') time, toteng, temp, press

!     write statistics file

        WRITE (9,'(1p,4e16.8)') time, toteng, temp, press

!     write trajectory file

        IF (ltraj .AND. MOD (istep, isample)==0) THEN
          
          WRITE (11,'(i10)') nbeads
          WRITE (11,'(f12.4,f12.7)') time, boxside

          IF (lpos) THEN

            DO i = 1, nbeads
              WRITE (11,'(a4,3f12.7)') 'Ar  ', x(i), y(i), z(i)
            END DO

          ELSE

            DO i = 1, nbeads
              WRITE (11,'(a4,3f12.7)') 'Ar  ', vx(i), vy(i), vz(i)
            END DO

          END IF

        END IF

      END DO

!     find and write calculation time

      CALL timecheck (timesim)
      WRITE (*,'("calculation time  = ",f12.3," seconds")') timesim
      WRITE (*,'("time per timestep = ",f12.6," seconds")') timesim/REAL(nstep, KIND=dp)

!     final average energy, temperature and pressure

      WRITE (*,*) 'final average energy, temperature and pressure:'
      WRITE (*,'(1x,i10,1p,3e16.8)') nav, totengav, tempav, pressav

!     write final configuration

      OPEN (8, file="XYZ")
      WRITE (8,'(i10,f12.7)') nbeads, boxside
      WRITE (8,'(a)') 'DPD Configuration'
      DO i = 1, nbeads
        WRITE (8,'(a4,3f12.7)') 'Ar  ', x(i), y(i), z(i)
      END DO

!     close all outputs and print end message

      IF (ltraj) CLOSE (11)      
      CLOSE (8)
      CLOSE (9)
      WRITE (*,*)'Job done'

END PROGRAM dpd

#ifdef FAST
      SUBROUTINE linkcell

!***********************************************************************
!
!     create link-cell lists for pairwise force calculations
!
!***********************************************************************

      USE variables
      IMPLICIT none

      INTEGER :: i, icell, ix, iy, iz, j

!     initialise link cell arrays

      lct = 0

!     calculate link cell indices

      DO i = 1, nbeads

        ix = INT ((x(i) + 0.5_dp * boxside) / wdthlnk)
        iy = INT ((y(i) + 0.5_dp * boxside) / wdthlnk)
        iz = INT ((z(i) + 0.5_dp * boxside) / wdthlnk)
        icell = 1 + ix + nlink * (iy + nlink * iz)
        j = lct (icell)
        lct (icell) = i
        link (i) = j

      END DO

      RETURN
      
      END SUBROUTINE linkcell

      SUBROUTINE dpd_forces

!***********************************************************************
!
!     calculate DPD conservative, dissipative and random forces using
!     link-cell lists to find particle pairs
!
!***********************************************************************

      USE variables
      IMPLICIT none

      INTEGER :: i, j, ix, iy, iz, jx, jy, jz, ibox, jbox, kc
      REAL(KIND=dp) :: dx, dy, dz, fcr, fdr, frr, pot, ftr, rrr, rsq, rdv

      INTEGER, DIMENSION(14) :: nix = (/ 0, 1, 1,-1, 0, 0,-1, 1,-1, 0, 1,-1, 0, 1 /)
      INTEGER, DIMENSION(14) :: niy = (/ 0, 0, 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1 /)
      INTEGER, DIMENSION(14) :: niz = (/ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 /)

      poteng = 0.0_dp
      virial = 0.0_dp

      fx = 0.0_dp
      fy = 0.0_dp
      fz = 0.0_dp

!     primary loops over all cells

      DO iz = 0, nlink-1

        DO iy = 0, nlink-1

          DO ix = 0, nlink-1

            ibox = 1 + ix + nlink * (iy + nlink * iz)

!     bypass cell if empty

            IF (lct (ibox)>0) THEN

!     secondary loop over neighbouring cells (including across periodic boundaries)

              DO kc = 1, 14

                jx = MOD(ix + nix (kc) + nlink, nlink)
                jy = MOD(iy + niy (kc) + nlink, nlink)
                jz = MOD(iz + niz (kc) + nlink, nlink)

!     index of neighbouring cell

                jbox = 1 + jx + nlink * (jy + nlink * jz)

!     bypass cell if empty

                IF (lct (jbox)>0) THEN

!     head of chain for ith subcell

                  i = lct (ibox)

!     loop over primary cell contents

                  DO WHILE (i/=0)

                    j = lct (jbox)

                    IF (ibox==jbox) j = link (i)
                    IF (j/=0) THEN

!     loop over secondary cell contents

                      DO WHILE (j/=0)

                        dx = x(j) - x(i)
                        dy = y(j) - y(i)
                        dz = z(j) - z(i)

!     calculate minimum image (across periodic boundaries)

                        dx = dx - boxside * NINT (dx/boxside)
                        dy = dy - boxside * NINT (dy/boxside)
                        dz = dz - boxside * NINT (dz/boxside)

!     calculate interparticle distance

                        rsq = dx * dx + dy * dy + dz * dz

!     spherical cutoff (rcut = 1.0_dp)

                        IF (rsq<1.0_dp) THEN

                          rrr = SQRT (rsq)

!     calculate dot product of vector and relative velocity

                          rdv = dx * (vx(j) - vx(i)) + dy * (vy(j) - vy(i)) + dz * (vz(j) - vz(i))

!     calculate forces for given particle pair

                          CALL calculate_dpd_forces (rrr, rsq, rdv, fcr, fdr, frr, pot)

!     accumulate potential and virial

                          poteng = poteng + pot
                          virial = virial - fcr * rsq

!     assign forces to particles

                          ftr = fcr + fdr + frr
                          fx(i) = fx(i) - ftr * dx
                          fy(i) = fy(i) - ftr * dy
                          fz(i) = fz(i) - ftr * dz
                          fx(j) = fx(j) + ftr * dx
                          fy(j) = fy(j) + ftr * dy
                          fz(j) = fz(j) + ftr * dz

                        END IF


                        j = link (j)

                      END DO

                    END IF

                    i = link (i)

                  END DO

                END IF

              END DO

            END IF

          END DO

        END DO

      END DO

      RETURN

      END SUBROUTINE dpd_forces
#else
      SUBROUTINE dpd_forces

!***********************************************************************
!     
!     calculate DPD conservative, dissipative and random forces
!     using brute force searching for particle pairs
!     
!***********************************************************************

      USE variables
      IMPLICIT none

      INTEGER :: i, j
      REAL(KIND=dp) :: dx, dy, dz, fcr, fdr, frr, pot, ftr, rrr, rsq, rdv

      poteng = 0.0_dp
      virial = 0.0_dp

      fx = 0.0_dp
      fy = 0.0_dp
      fz = 0.0_dp

!     loop over all particle pairs
 
      DO i = 1, nbeads-1

        DO j = i+1, nbeads

          dx = x(j) - x(i)
          dy = y(j) - y(i)
          dz = z(j) - z(i)

!     calculate minimum image (across periodic boundaries)

          dx = dx - boxside * NINT (dx/boxside)
          dy = dy - boxside * NINT (dy/boxside)
          dz = dz - boxside * NINT (dz/boxside)

!     calculate interparticle distance

          rsq = dx * dx + dy * dy + dz * dz

!     spherical cutoff (rcut = 1.0_dp)

          IF (rsq<1.0_dp) THEN

            rrr = SQRT (rsq)

!     calculate dot product of vector and relative velocity

            rdv = dx * (vx(j) - vx(i)) + dy * (vy(j) - vy(i)) + dz * (vz(j) - vz(i))

!     calculate forces for given particle pair

            CALL calculate_dpd_forces (rrr, rsq, rdv, fcr, fdr, frr, pot)

!     accumulate potential and virial

            poteng = poteng + pot
            virial = virial - fcr * rsq

!     assign forces to particles

            ftr = fcr + fdr + frr
            fx(i) = fx(i) - ftr * dx
            fy(i) = fy(i) - ftr * dy
            fz(i) = fz(i) - ftr * dz
            fx(j) = fx(j) + ftr * dx
            fy(j) = fy(j) + ftr * dy
            fz(j) = fz(j) + ftr * dz

          END IF

        END DO

      END DO

      RETURN

      END SUBROUTINE dpd_forces
#endif

      SUBROUTINE calculate_dpd_forces (rrr, rsq, rdv, fcr, fdr, frr, pot)

!***********************************************************************
!
!     calculate DPD conservative, dissipative and random forces
!     for a given particle pair (inputs are distance, squared distance,
!     and dot product of vector and relative velocity, outputs are
!     forces and interaction potential)
!
!***********************************************************************

      USE variables
      IMPLICIT none

      REAL(KIND=dp), INTENT (IN) :: rrr, rsq, rdv
      REAL(KIND=dp), INTENT (OUT) :: fcr, fdr, frr, pot

      REAL(KIND=dp) :: wrr, zeta

      fcr = 0.0_dp
      fdr = 0.0_dp
      frr = 0.0_dp
      pot = 0.0_dp

      IF (rrr>1.0e-10_dp) THEN

!     screening function (for conservative force)

         wrr = MAX ((bij - rrr), 0.0_dp)

!     potential calculation

         pot = 0.5_dp * aij * wrr * wrr

!     conservative interaction force

         fcr = aij * wrr / rrr

!     dissipative force (rdv = dot product of vector and relative velocity)

         fdr = -gamma * rdv * wrr * wrr / rsq

!     random force (zeta = gaussian random number)

         CALL RANDOM_NUMBER (zeta)
         zeta = fac12 * (zeta - 0.5_dp)
         frr = zeta * sigma * wrr / rrr

      END IF

      RETURN
      END SUBROUTINE calculate_dpd_forces

      SUBROUTINE timecheck (time)

!***********************************************************************
!
!     timing routine (gives time elapsed in seconds)
!
!***********************************************************************

      USE variables
      IMPLICIT none

      INTEGER :: ticks
      INTEGER, SAVE :: init = 0, tickrate, maxtick
      REAL(KIND=dp) :: time
      REAL(KIND=dp), SAVE :: tzero

      IF (init==0) THEN

        init = 1
        CALL SYSTEM_CLOCK (count_max=maxtick,count_rate=tickrate)
        CALL SYSTEM_CLOCK (ticks)
        ticks = MOD (ticks+maxtick, maxtick)
        tzero = REAL (ticks, KIND=dp)/ REAL (tickrate, KIND=dp)

      END IF

      CALL SYSTEM_CLOCK (ticks)
      ticks = MOD (ticks+maxtick, maxtick)
      time = REAL (ticks, KIND=dp)/ REAL (tickrate, KIND=dp) - tzero

      RETURN
      END SUBROUTINE timecheck
